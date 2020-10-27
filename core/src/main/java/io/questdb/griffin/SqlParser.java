/*******************************************************************************
 *     ___                  _   ____  ____
 *    / _ \ _   _  ___  ___| |_|  _ \| __ )
 *   | | | | | | |/ _ \/ __| __| | | |  _ \
 *   | |_| | |_| |  __/\__ \ |_| |_| | |_) |
 *    \__\_\\__,_|\___||___/\__|____/|____/
 *
 *  Copyright (c) 2014-2019 Appsicle
 *  Copyright (c) 2019-2020 QuestDB
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 ******************************************************************************/

package io.questdb.griffin;

import io.questdb.cairo.CairoConfiguration;
import io.questdb.cairo.ColumnType;
import io.questdb.cairo.PartitionBy;
import io.questdb.cairo.TableUtils;
import io.questdb.cairo.sql.RecordMetadata;
import io.questdb.griffin.model.AnalyticColumn;
import io.questdb.griffin.model.ColumnCastModel;
import io.questdb.griffin.model.CopyModel;
import io.questdb.griffin.model.CreateTableModel;
import io.questdb.griffin.model.ExecutionModel;
import io.questdb.griffin.model.ExpressionNode;
import io.questdb.griffin.model.InsertModel;
import io.questdb.griffin.model.QueryColumn;
import io.questdb.griffin.model.QueryModel;
import io.questdb.griffin.model.RenameTableModel;
import io.questdb.griffin.model.UpdateModel;
import io.questdb.griffin.model.WithClauseModel;
import io.questdb.std.CharSequenceHashSet;
import io.questdb.std.CharSequenceObjHashMap;
import io.questdb.std.Chars;
import io.questdb.std.GenericLexer;
import io.questdb.std.LowerCaseAsciiCharSequenceHashSet;
import io.questdb.std.LowerCaseAsciiCharSequenceIntHashMap;
import io.questdb.std.Numbers;
import io.questdb.std.NumericException;
import io.questdb.std.ObjList;
import io.questdb.std.ObjectPool;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import static io.questdb.griffin.SqlKeywords.isAllKeyword;
import static io.questdb.griffin.SqlKeywords.isAsKeyword;
import static io.questdb.griffin.SqlKeywords.isAscKeyword;
import static io.questdb.griffin.SqlKeywords.isByKeyword;
import static io.questdb.griffin.SqlKeywords.isCacheKeyword;
import static io.questdb.griffin.SqlKeywords.isCapacityKeyword;
import static io.questdb.griffin.SqlKeywords.isCaseKeyword;
import static io.questdb.griffin.SqlKeywords.isCastKeyword;
import static io.questdb.griffin.SqlKeywords.isCopyKeyword;
import static io.questdb.griffin.SqlKeywords.isCountKeyword;
import static io.questdb.griffin.SqlKeywords.isCreateKeyword;
import static io.questdb.griffin.SqlKeywords.isDescKeyword;
import static io.questdb.griffin.SqlKeywords.isDistinctKeyword;
import static io.questdb.griffin.SqlKeywords.isExceptKeyword;
import static io.questdb.griffin.SqlKeywords.isFillKeyword;
import static io.questdb.griffin.SqlKeywords.isFromKeyword;
import static io.questdb.griffin.SqlKeywords.isGroupKeyword;
import static io.questdb.griffin.SqlKeywords.isHeaderKeyword;
import static io.questdb.griffin.SqlKeywords.isIndexKeyword;
import static io.questdb.griffin.SqlKeywords.isInsertKeyword;
import static io.questdb.griffin.SqlKeywords.isIntersectKeyword;
import static io.questdb.griffin.SqlKeywords.isJoinKeyword;
import static io.questdb.griffin.SqlKeywords.isLatestKeyword;
import static io.questdb.griffin.SqlKeywords.isLimitKeyword;
import static io.questdb.griffin.SqlKeywords.isNoCacheKeyword;
import static io.questdb.griffin.SqlKeywords.isOnKeyword;
import static io.questdb.griffin.SqlKeywords.isOrderKeyword;
import static io.questdb.griffin.SqlKeywords.isOverKeyword;
import static io.questdb.griffin.SqlKeywords.isPartitionKeyword;
import static io.questdb.griffin.SqlKeywords.isRenameKeyword;
import static io.questdb.griffin.SqlKeywords.isSampleKeyword;
import static io.questdb.griffin.SqlKeywords.isSelectKeyword;
import static io.questdb.griffin.SqlKeywords.isTimestampKeyword;
import static io.questdb.griffin.SqlKeywords.isTrueKeyword;
import static io.questdb.griffin.SqlKeywords.isUnionKeyword;
import static io.questdb.griffin.SqlKeywords.isUpdateKeyword;
import static io.questdb.griffin.SqlKeywords.isValuesKeyword;
import static io.questdb.griffin.SqlKeywords.isWhereKeyword;
import static io.questdb.griffin.SqlKeywords.isWithKeyword;

public final class SqlParser {

    public static final int MAX_ORDER_BY_COLUMNS = 1560;
    private static final LowerCaseAsciiCharSequenceHashSet tableAliasStop = new LowerCaseAsciiCharSequenceHashSet();
    private static final LowerCaseAsciiCharSequenceHashSet columnAliasStop = new LowerCaseAsciiCharSequenceHashSet();
    private static final LowerCaseAsciiCharSequenceHashSet groupByStopSet = new LowerCaseAsciiCharSequenceHashSet();
    private static final LowerCaseAsciiCharSequenceIntHashMap joinStartSet = new LowerCaseAsciiCharSequenceIntHashMap();
    private static final LowerCaseAsciiCharSequenceHashSet setOperations = new LowerCaseAsciiCharSequenceHashSet();

    static {
        tableAliasStop.add("where");
        tableAliasStop.add("latest");
        tableAliasStop.add("join");
        tableAliasStop.add("inner");
        tableAliasStop.add("left");
        tableAliasStop.add("outer");
        tableAliasStop.add("asof");
        tableAliasStop.add("splice");
        tableAliasStop.add("lt");
        tableAliasStop.add("cross");
        tableAliasStop.add("sample");
        tableAliasStop.add("order");
        tableAliasStop.add("on");
        tableAliasStop.add("timestamp");
        tableAliasStop.add("limit");
        tableAliasStop.add(")");
        tableAliasStop.add(";");
        tableAliasStop.add("union");
        tableAliasStop.add("group");
        tableAliasStop.add("except");
        tableAliasStop.add("intersect");
        //
        columnAliasStop.add("from");
        columnAliasStop.add(",");
        columnAliasStop.add("over");
        //
        groupByStopSet.add("order");
        groupByStopSet.add(")");
        groupByStopSet.add(",");

        joinStartSet.put("left", QueryModel.JOIN_INNER);
        joinStartSet.put("join", QueryModel.JOIN_INNER);
        joinStartSet.put("inner", QueryModel.JOIN_INNER);
        joinStartSet.put("outer", QueryModel.JOIN_OUTER);
        joinStartSet.put("cross", QueryModel.JOIN_CROSS);
        joinStartSet.put("asof", QueryModel.JOIN_ASOF);
        joinStartSet.put("splice", QueryModel.JOIN_SPLICE);
        joinStartSet.put("lt", QueryModel.JOIN_LT);
        joinStartSet.put(",", QueryModel.JOIN_CROSS);
        //
        setOperations.add("union");
        setOperations.add("except");
        setOperations.add("intersect");
    }


    private final ObjectPool<ExpressionNode> expressionNodePool;
    private final ExpressionTreeBuilder expressionTreeBuilder;
    private final ObjectPool<QueryModel> queryModelPool;
    private final ObjectPool<QueryColumn> queryColumnPool;
    private final ObjectPool<AnalyticColumn> analyticColumnPool;
    private final ObjectPool<CreateTableModel> createTableModelPool;
    private final ObjectPool<ColumnCastModel> columnCastModelPool;
    private final ObjectPool<RenameTableModel> renameTableModelPool;
    private final ObjectPool<WithClauseModel> withClauseModelPool;
    private final ObjectPool<InsertModel> insertModelPool;
    private final ObjectPool<CopyModel> copyModelPool;
    private final ObjectPool<UpdateModel> updateModelObjectPool;
    private final ExpressionParser expressionParser;
    private final CairoConfiguration configuration;
    private final PostOrderTreeTraversalAlgo traversalAlgo;
    private final ObjList<ExpressionNode> tempExprNodes = new ObjList<>();
    private final CharacterStore characterStore;
    private final SqlOptimiser optimiser;
    private final PostOrderTreeTraversalAlgo.Visitor rewriteCase0Ref = this::rewriteCase0;
    private final PostOrderTreeTraversalAlgo.Visitor rewriteCount0Ref = this::rewriteCount0;
    private boolean subQueryMode = false;

    SqlParser(
            CairoConfiguration configuration,
            SqlOptimiser optimiser,
            CharacterStore characterStore,
            ObjectPool<ExpressionNode> expressionNodePool,
            ObjectPool<QueryColumn> queryColumnPool,
            ObjectPool<QueryModel> queryModelPool,
            PostOrderTreeTraversalAlgo traversalAlgo
    ) {
        this.expressionNodePool = expressionNodePool;
        this.queryModelPool = queryModelPool;
        this.queryColumnPool = queryColumnPool;
        this.expressionTreeBuilder = new ExpressionTreeBuilder();
        this.analyticColumnPool = new ObjectPool<>(AnalyticColumn.FACTORY, configuration.getAnalyticColumnPoolCapacity());
        this.createTableModelPool = new ObjectPool<>(CreateTableModel.FACTORY, configuration.getCreateTableModelPoolCapacity());
        this.columnCastModelPool = new ObjectPool<>(ColumnCastModel.FACTORY, configuration.getColumnCastModelPoolCapacity());
        this.renameTableModelPool = new ObjectPool<>(RenameTableModel.FACTORY, configuration.getRenameTableModelPoolCapacity());
        this.withClauseModelPool = new ObjectPool<>(WithClauseModel.FACTORY, configuration.getWithClauseModelPoolCapacity());
        this.insertModelPool = new ObjectPool<>(InsertModel.FACTORY, configuration.getInsertPoolCapacity());
        this.copyModelPool = new ObjectPool<>(CopyModel.FACTORY, configuration.getCopyPoolCapacity());
        this.updateModelObjectPool = new ObjectPool<>(UpdateModel.FACTORY, configuration.getUpdatePoolCapacity());
        this.configuration = configuration;
        this.traversalAlgo = traversalAlgo;
        this.characterStore = characterStore;
        this.optimiser = optimiser;
        this.expressionParser = new ExpressionParser(expressionNodePool, this, characterStore);
    }

    private static SqlException err(GenericLexer lexer, String msg) {
        return SqlException.$(lexer.lastTokenPosition(), msg);
    }

    private static SqlException errUnexpected(GenericLexer lexer, CharSequence token) {
        return SqlException.unexpectedToken(lexer.lastTokenPosition(), token);
    }

    private void assertNotDot(GenericLexer lexer, CharSequence tok) throws SqlException {
        if (Chars.indexOf(tok, '.') != -1) {
            throw SqlException.$(lexer.lastTokenPosition(), "'.' is not allowed here");
        }
    }

    void clear() {
        queryModelPool.clear();
        queryColumnPool.clear();
        expressionNodePool.clear();
        analyticColumnPool.clear();
        createTableModelPool.clear();
        columnCastModelPool.clear();
        renameTableModelPool.clear();
        withClauseModelPool.clear();
        subQueryMode = false;
        characterStore.clear();
        insertModelPool.clear();
        expressionTreeBuilder.reset();
        copyModelPool.clear();
    }

    private CharSequence createColumnAlias(ExpressionNode node, QueryModel model) {
        return SqlUtil.createColumnAlias(
                characterStore,
                GenericLexer.unquote(node.token),
                Chars.indexOf(node.token, '.'),
                model.getAliasToColumnMap()
        );
    }

    private ExpressionNode expectExpr(GenericLexer lexer) throws SqlException {
        ExpressionNode n = expr(lexer, (QueryModel) null);
        if (n == null) {
            throw SqlException.$(lexer.getUnparsed() == null ? lexer.getPosition() : lexer.lastTokenPosition(), "Expression expected");
        }
        return n;
    }

    private int expectInt(GenericLexer lexer) throws SqlException {
        CharSequence tok = tok(lexer, "integer");
        boolean negative;
        if (Chars.equals(tok, '-')) {
            negative = true;
            tok = tok(lexer, "integer");
        } else {
            negative = false;
        }
        try {
            int result = Numbers.parseInt(tok);
            return negative ? -result : result;
        } catch (NumericException e) {
            throw err(lexer, "bad integer");
        }
    }

    private ExpressionNode expectLiteral(GenericLexer lexer) throws SqlException {
        CharSequence tok = tok(lexer, "literal");
        int pos = lexer.lastTokenPosition();
        validateLiteral(pos, tok);
        return nextLiteral(GenericLexer.immutableOf(GenericLexer.unquote(tok)), pos);
    }

    private CharSequence expectTableNameOrSubQuery(GenericLexer lexer) throws SqlException {
        return tok(lexer, "table name or sub-query");
    }

    private void expectTok(GenericLexer lexer, CharSequence tok, CharSequence expected) throws SqlException {
        if (tok == null || !Chars.equalsLowerCaseAscii(tok, expected)) {
            throw SqlException.position(lexer.lastTokenPosition()).put('\'').put(expected).put("' expected");
        }
    }

    private void expectTok(GenericLexer lexer, CharSequence expected) throws SqlException {
        CharSequence tok = optTok(lexer);
        if (tok == null) {
            throw SqlException.position(lexer.getPosition()).put('\'').put(expected).put("' expected");
        }
        expectTok(lexer, tok, expected);
    }

    private void expectTok(GenericLexer lexer, char expected) throws SqlException {
        CharSequence tok = optTok(lexer);
        if (tok == null) {
            throw SqlException.position(lexer.getPosition()).put('\'').put(expected).put("' expected");
        }
        expectTok(tok, lexer.lastTokenPosition(), expected);
    }

    private void expectTok(CharSequence tok, int pos, char expected) throws SqlException {
        if (tok == null || !Chars.equals(tok, expected)) {
            throw SqlException.position(pos).put('\'').put(expected).put("' expected");
        }
    }

    ExpressionNode expr(GenericLexer lexer, QueryModel model) throws SqlException {
        try {
            expressionTreeBuilder.pushModel(model);
            expressionParser.parseExpr(lexer, expressionTreeBuilder);
            return rewriteCase(rewriteCount(expressionTreeBuilder.poll()));
        } catch (SqlException e) {
            expressionTreeBuilder.reset();
            throw e;
        } finally {
            expressionTreeBuilder.popModel();
        }
    }

    // test only
    void expr(GenericLexer lexer, ExpressionParserListener listener) throws SqlException {
        expressionParser.parseExpr(lexer, listener);
    }

    private int getCreateTableColumnIndex(CreateTableModel model, CharSequence columnName, int position) throws SqlException {
        int index = model.getColumnIndex(columnName);
        if (index == -1) {
            throw SqlException.invalidColumn(position, columnName);
        }
        return index;
    }

    private boolean isFieldTerm(CharSequence tok) {
        return Chars.equals(tok, ')') || Chars.equals(tok, ',');
    }

    private ExpressionNode literal(GenericLexer lexer, CharSequence name) {
        return literal(name, lexer.lastTokenPosition());
    }

    private ExpressionNode literal(CharSequence name, int position) {
        // this can never be null in its current contexts
        // every time this function is called is after lexer.unparse(), which ensures non-null token.
        return expressionNodePool.next().of(ExpressionNode.LITERAL, GenericLexer.unquote(name), 0, position);
    }

    private ExpressionNode nextLiteral(CharSequence token, int position) {
        return SqlUtil.nextLiteral(expressionNodePool, token, position);
    }

    private CharSequence notTermTok(GenericLexer lexer) throws SqlException {
        CharSequence tok = tok(lexer, "')' or ','");
        if (isFieldTerm(tok)) {
            throw err(lexer, "missing column definition");
        }
        return tok;
    }

    private CharSequence optTok(GenericLexer lexer) {
        CharSequence tok = SqlUtil.fetchNext(lexer);
        if (tok == null || (subQueryMode && Chars.equals(tok, ')'))) {
            return null;
        }
        return tok;
    }

    ExecutionModel parse(GenericLexer lexer, SqlExecutionContext executionContext) throws SqlException {
        CharSequence tok = tok(lexer, "'create', 'rename', 'select' or 'update'");

        if (isSelectKeyword(tok)) {
            return parseSelect(lexer);
        }

        if (isCreateKeyword(tok)) {
            return parseCreateStatement(lexer, executionContext);
        }

        if (isRenameKeyword(tok)) {
            return parseRenameStatement(lexer);
        }

        if (isInsertKeyword(tok)) {
            return parseInsert(lexer);
        }

        if (isCopyKeyword(tok)) {
            return parseCopy(lexer);
        }

        if(isUpdateKeyword(tok)) {
            return parseUpdate(lexer);
        }

        return parseSelect(lexer);
    }

    QueryModel parseAsSubQuery(GenericLexer lexer, @Nullable CharSequenceObjHashMap<WithClauseModel> withClauses) throws SqlException {
        QueryModel model;
        this.subQueryMode = true;
        try {
            model = parseDml(lexer, withClauses);
        } finally {
            this.subQueryMode = false;
        }
        return model;
    }

    private QueryModel parseAsSubQueryAndExpectClosingBrace(GenericLexer lexer, CharSequenceObjHashMap<WithClauseModel> withClauses) throws SqlException {
        final QueryModel model = parseAsSubQuery(lexer, withClauses);
        expectTok(lexer, ')');
        return model;
    }

    private ExecutionModel parseCopy(GenericLexer lexer) throws SqlException {
        if (configuration.getInputRoot() == null) {
            throw SqlException.$(lexer.lastTokenPosition(), "COPY is disabled ['cairo.sql.copy.root' is not set?]");
        }
        ExpressionNode tableName = expectExpr(lexer);
        CharSequence tok = tok(lexer, "'from' or 'to'");

        if (isFromKeyword(tok)) {
            final ExpressionNode fileName = expectExpr(lexer);
            if (fileName.token.length() < 3 && Chars.startsWith(fileName.token, '\'')) {
                throw SqlException.$(fileName.position, "file name expected");
            }
            CopyModel model = copyModelPool.next();
            model.setTableName(tableName);
            model.setFileName(fileName);

            tok = optTok(lexer);
            if (tok != null && isWithKeyword(tok)) {
                tok = tok(lexer, "copy option");
                while (tok != null) {
                    if (isHeaderKeyword(tok)) {
                        model.setHeader(isTrueKeyword(tok(lexer, "'true' or 'false'")));
                        tok = optTok(lexer);
                    } else {
                        throw SqlException.$(lexer.lastTokenPosition(), "unexpected option");
                    }
                }
            }
            return model;
        }
        throw SqlException.$(lexer.lastTokenPosition(), "'from' expected");
    }

    private ExecutionModel parseCreateStatement(GenericLexer lexer, SqlExecutionContext executionContext) throws SqlException {
        expectTok(lexer, "table");
        return parseCreateTable(lexer, executionContext);
    }

    private ExecutionModel parseCreateTable(GenericLexer lexer, SqlExecutionContext executionContext) throws SqlException {
        final CreateTableModel model = createTableModelPool.next();
        final CharSequence tableName = tok(lexer, "table name");
        model.setName(nextLiteral(GenericLexer.assertNoDotsAndSlashes(GenericLexer.unquote(tableName), lexer.lastTokenPosition()), lexer.lastTokenPosition()));

        CharSequence tok = tok(lexer, "'(' or 'as'");

        if (Chars.equals(tok, '(')) {
            lexer.unparse();
            parseCreateTableColumns(lexer, model);
        } else if (isAsKeyword(tok)) {
            parseCreateTableAsSelect(lexer, model, executionContext);
        } else {
            throw errUnexpected(lexer, tok);
        }

        while ((tok = optTok(lexer)) != null && Chars.equals(tok, ',')) {
            tok = tok(lexer, "'index' or 'cast'");
            if (isIndexKeyword(tok)) {
                parseCreateTableIndexDef(lexer, model);
            } else if (isCastKeyword(tok)) {
                parseCreateTableCastDef(lexer, model);
            } else {
                throw errUnexpected(lexer, tok);
            }
        }

        ExpressionNode timestamp = parseTimestamp(lexer, tok);
        if (timestamp != null) {
            // ignore index, validate column
            getCreateTableColumnIndex(model, timestamp.token, timestamp.position);
            model.setTimestamp(timestamp);
            tok = optTok(lexer);
        }

        ExpressionNode partitionBy = parseCreateTablePartition(lexer, tok);
        if (partitionBy != null) {
            if (PartitionBy.fromString(partitionBy.token) == -1) {
                throw SqlException.$(partitionBy.position, "'NONE', 'DAY', 'MONTH' or 'YEAR' expected");
            }
            model.setPartitionBy(partitionBy);
            tok = optTok(lexer);
        }

        if (tok == null || Chars.equals(tok, ';')) {
            return model;
        }
        throw errUnexpected(lexer, tok);
    }

    private void parseCreateTableAsSelect(GenericLexer lexer, CreateTableModel model, SqlExecutionContext executionContext) throws SqlException {
        expectTok(lexer, '(');
        QueryModel queryModel = optimiser.optimise(parseDml(lexer, null), executionContext);
        ObjList<QueryColumn> columns = queryModel.getBottomUpColumns();
        assert columns.size() > 0;

        // we do not know types of columns at this stage
        // compiler must put table together using query metadata.
        for (int i = 0, n = columns.size(); i < n; i++) {
            model.addColumn(columns.getQuick(i).getName(), -1, configuration.getDefaultSymbolCapacity());
        }

        model.setQueryModel(queryModel);
        expectTok(lexer, ')');
    }

    private void parseCreateTableCastDef(GenericLexer lexer, CreateTableModel model) throws SqlException {
        if (model.getQueryModel() == null) {
            throw SqlException.$(lexer.lastTokenPosition(), "cast is only supported in 'create table as ...' context");
        }
        expectTok(lexer, '(');
        ColumnCastModel columnCastModel = columnCastModelPool.next();

        final ExpressionNode columnName = expectLiteral(lexer);
        columnCastModel.setName(columnName);
        expectTok(lexer, "as");

        final ExpressionNode columnType = expectLiteral(lexer);
        final int type = toColumnType(lexer, columnType.token);
        columnCastModel.setType(type, columnName.position, columnType.position);

        if (type == ColumnType.SYMBOL) {
            CharSequence tok = tok(lexer, "'capacity', 'nocache', 'cache', 'index' or ')'");

            int symbolCapacity;
            int capacityPosition;
            if (isCapacityKeyword(tok)) {
                capacityPosition = lexer.getPosition();
                columnCastModel.setSymbolCapacity(symbolCapacity = parseSymbolCapacity(lexer));
                tok = tok(lexer, "'nocache', 'cache', 'index' or ')'");
            } else {
                columnCastModel.setSymbolCapacity(configuration.getDefaultSymbolCapacity());
                symbolCapacity = -1;
                capacityPosition = -1;
            }

            final boolean cached;
            if (isNoCacheKeyword(tok)) {
                cached = false;
            } else if (isCacheKeyword(tok)) {
                cached = true;
            } else {
                cached = configuration.getDefaultSymbolCacheFlag();
                lexer.unparse();
            }

            columnCastModel.setSymbolCacheFlag(cached);

            if (cached && symbolCapacity != -1) {
                assert capacityPosition != -1;
                TableUtils.validateSymbolCapacityCached(true, symbolCapacity, capacityPosition);
            }

            // parse index clause

            tok = tok(lexer, "')', or 'index'");

            if (isIndexKeyword(tok)) {
                columnCastModel.setIndexed(true);
                tok = tok(lexer, "')', or 'capacity'");

                if (isCapacityKeyword(tok)) {
                    int errorPosition = lexer.getPosition();
                    int indexValueBlockSize = expectInt(lexer);
                    TableUtils.validateIndexValueBlockSize(errorPosition, indexValueBlockSize);
                    columnCastModel.setIndexValueBlockSize(Numbers.ceilPow2(indexValueBlockSize));
                } else {
                    columnCastModel.setIndexValueBlockSize(configuration.getIndexValueBlockSize());
                }
            } else {
                columnCastModel.setIndexed(false);
                lexer.unparse();
            }
        }

        expectTok(lexer, ')');

        if (!model.addColumnCastModel(columnCastModel)) {
            throw SqlException.$(columnCastModel.getName().position, "duplicate cast");
        }
    }

    private void parseCreateTableColumns(GenericLexer lexer, CreateTableModel model) throws SqlException {
        expectTok(lexer, '(');

        while (true) {
            final int position = lexer.lastTokenPosition();
            final CharSequence name = GenericLexer.immutableOf(GenericLexer.unquote(notTermTok(lexer)));
            final int type = toColumnType(lexer, notTermTok(lexer));

            if (SqlKeywords.isInvalidColumnName(name)) {
                throw SqlException.$(position, " new column name contains invalid characters");
            }

            if (!model.addColumn(name, type, configuration.getDefaultSymbolCapacity())) {
                throw SqlException.$(position, "Duplicate column");
            }

            CharSequence tok;
            if (type == ColumnType.SYMBOL) {
                tok = tok(lexer, "'capacity', 'nocache', 'cache', 'index' or ')'");

                int symbolCapacity;
                if (isCapacityKeyword(tok)) {
                    // when capacity is not set explicitly it will default via configuration
                    model.symbolCapacity(symbolCapacity = parseSymbolCapacity(lexer));
                    tok = tok(lexer, "'nocache', 'cache', 'index' or ')'");
                } else {
                    symbolCapacity = -1;
                }

                final boolean cached;
                if (isNoCacheKeyword(tok)) {
                    cached = false;
                } else if (isCacheKeyword(tok)) {
                    cached = true;
                } else {
                    cached = configuration.getDefaultSymbolCacheFlag();
                    lexer.unparse();
                }
                model.cached(cached);
                if (cached && symbolCapacity != -1) {
                    TableUtils.validateSymbolCapacityCached(true, symbolCapacity, lexer.lastTokenPosition());
                }
                tok = parseCreateTableInlineIndexDef(lexer, model);
            } else {
                tok = null;
            }

            if (tok == null) {
                tok = tok(lexer, "',' or ')'");
            }

            if (Chars.equals(tok, ')')) {
                break;
            }

            if (!Chars.equals(tok, ',')) {
                throw err(lexer, "',' or ')' expected");
            }
        }
    }

    private void parseCreateTableIndexDef(GenericLexer lexer, CreateTableModel model) throws SqlException {
        expectTok(lexer, '(');
        final int columnIndex = getCreateTableColumnIndex(model, expectLiteral(lexer).token, lexer.lastTokenPosition());

        if (isCapacityKeyword(tok(lexer, "'capacity'"))) {
            int errorPosition = lexer.getPosition();
            int indexValueBlockSize = expectInt(lexer);
            TableUtils.validateIndexValueBlockSize(errorPosition, indexValueBlockSize);
            model.setIndexFlags(columnIndex, true, Numbers.ceilPow2(indexValueBlockSize));
        } else {
            model.setIndexFlags(columnIndex, true, configuration.getIndexValueBlockSize());
            lexer.unparse();
        }
        expectTok(lexer, ')');
    }

    private CharSequence parseCreateTableInlineIndexDef(GenericLexer lexer, CreateTableModel model) throws SqlException {
        CharSequence tok = tok(lexer, "')', or 'index'");

        if (isFieldTerm(tok)) {
            model.setIndexFlags(false, configuration.getIndexValueBlockSize());
            return tok;
        }

        expectTok(lexer, tok, "index");

        if (isFieldTerm(tok = tok(lexer, ") | , expected"))) {
            model.setIndexFlags(true, configuration.getIndexValueBlockSize());
            return tok;
        }

        expectTok(lexer, tok, "capacity");

        int errorPosition = lexer.getPosition();
        int indexValueBlockSize = expectInt(lexer);
        TableUtils.validateIndexValueBlockSize(errorPosition, indexValueBlockSize);
        model.setIndexFlags(true, Numbers.ceilPow2(indexValueBlockSize));
        return null;
    }

    private void expectBy(GenericLexer lexer) throws SqlException {
        CharSequence tok = optTok(lexer);
        if (tok == null || !isByKeyword(tok)) {
            throw SqlException.$((lexer.getPosition()), "'by' expected");
        }
    }

    private ExpressionNode parseCreateTablePartition(GenericLexer lexer, CharSequence tok) throws SqlException {
        if (tok != null && isPartitionKeyword(tok)) {
            expectTok(lexer, "by");
            return expectLiteral(lexer);
        }
        return null;
    }

    private QueryModel parseDml(GenericLexer lexer, @Nullable CharSequenceObjHashMap<WithClauseModel> withClauses) throws SqlException {
        QueryModel model = null;
        QueryModel prevModel = null;
        while (true) {

            QueryModel unionModel = parseDml0(lexer, prevModel != null ? prevModel.getWithClauses() : withClauses);
            if (prevModel == null) {
                model = unionModel;
                prevModel = model;
            } else {
                prevModel.setUnionModel(unionModel);
                prevModel = unionModel;
            }

            CharSequence tok = optTok(lexer);
            if (tok == null || setOperations.excludes(tok)) {
                lexer.unparse();
                return model;
            }

            if (isUnionKeyword(tok)) {
                tok = tok(lexer, "all or select");
                if (isAllKeyword(tok)) {
                    if (!model.isDistinct()) {
                        prevModel.setSetOperationType(QueryModel.SET_OPERATION_UNION_ALL);
                    } else {
                        prevModel.setSetOperationType(QueryModel.SET_OPERATION_UNION);
                    }
                } else {
                    prevModel.setSetOperationType(QueryModel.SET_OPERATION_UNION);
                    lexer.unparse();
                }
            }

            if (isExceptKeyword(tok)) {
                prevModel.setSetOperationType(QueryModel.SET_OPERATION_EXCEPT);
            }

            if (isIntersectKeyword(tok)) {
                prevModel.setSetOperationType(QueryModel.SET_OPERATION_INTERSECT);
            }
        }
    }

    private void parseSelectClause(GenericLexer lexer, QueryModel model) throws SqlException {
        CharSequence tok = tok(lexer, "[distinct] column");

        ExpressionNode expr;
        if (isDistinctKeyword(tok)) {
            model.setDistinct(true);
        } else {
            lexer.unparse();
        }
        while (true) {

            tok = tok(lexer, "column");
            if (Chars.equals(tok, '*')) {
                expr = nextLiteral(GenericLexer.immutableOf(tok), lexer.lastTokenPosition());
            } else {
                // cut off some obvious errors
                if (isFromKeyword(tok)) {
                    throw SqlException.$(lexer.getPosition(), "column name expected");
                }

                if (isSelectKeyword(tok)) {
                    throw SqlException.$(lexer.getPosition(), "reserved name");
                }

                lexer.unparse();
                expr = expr(lexer, model);

                if (expr == null) {
                    throw SqlException.$(lexer.lastTokenPosition(), "missing expression");
                }

                if (Chars.endsWith(expr.token, '.')) {
                    throw SqlException.$(expr.position + expr.token.length(), "'*' or column name expected");
                }
            }

            final CharSequence alias;

            tok = optTok(lexer);

            if (tok != null && Chars.equals(tok, ';')) {
                alias = createColumnAlias(expr, model);
                tok = optTok(lexer);
            } else if (tok != null && columnAliasStop.excludes(tok)) {
                assertNotDot(lexer, tok);

                if (isAsKeyword(tok)) {
                    alias = GenericLexer.unquote(GenericLexer.immutableOf(tok(lexer, "alias")));
                } else {
                    alias = GenericLexer.immutableOf(tok);
                }
                tok = optTok(lexer);
            } else {
                alias = createColumnAlias(expr, model);
            }

            if (tok != null && isOverKeyword(tok)) {
                // analytic
                expectTok(lexer, '(');

                AnalyticColumn col = analyticColumnPool.next().of(alias, expr);
                tok = tok(lexer, "'");

                if (isPartitionKeyword(tok)) {
                    expectTok(lexer, "by");

                    ObjList<ExpressionNode> partitionBy = col.getPartitionBy();

                    do {
                        partitionBy.add(expectLiteral(lexer));
                        tok = tok(lexer, "'order' or ')'");
                    } while (Chars.equals(tok, ','));
                }

                if (isOrderKeyword(tok)) {
                    expectTok(lexer, "by");

                    do {
                        ExpressionNode e = expectLiteral(lexer);
                        tok = tok(lexer, "'asc' or 'desc'");

                        if (isDescKeyword(tok)) {
                            col.addOrderBy(e, QueryModel.ORDER_DIRECTION_DESCENDING);
                            tok = tok(lexer, "',' or ')'");
                        } else {
                            col.addOrderBy(e, QueryModel.ORDER_DIRECTION_ASCENDING);
                            if (isAscKeyword(tok)) {
                                tok = tok(lexer, "',' or ')'");
                            }
                        }
                    } while (Chars.equals(tok, ','));
                }
                expectTok(tok, lexer.lastTokenPosition(), ')');
                model.addBottomUpColumn(col);
                tok = tok(lexer, "'from' or ','");
            } else {
                if (expr.type == ExpressionNode.QUERY) {
                    throw SqlException.$(expr.position, "query is not expected, did you mean column?");
                }
                model.addBottomUpColumn(queryColumnPool.next().of(alias, expr));
            }

            if (tok == null) {
                lexer.unparse();
                break;
            }

            if (isFromKeyword(tok)) {
                lexer.unparse();
                break;
            }

            if (setOperations.contains(tok)) {
                lexer.unparse();
                break;
            }

            if (!Chars.equals(tok, ',')) {
                throw err(lexer, "',' or 'from' expected");
            }
        }
    }

    @NotNull
    private QueryModel parseDml0(GenericLexer lexer, @Nullable CharSequenceObjHashMap<WithClauseModel> parentWithClauses) throws SqlException {
        CharSequence tok;
        final int modelPosition = lexer.getPosition();

        QueryModel model = queryModelPool.next();
        model.setModelPosition(modelPosition);
        if (parentWithClauses != null) {
            model.addWithClauses(parentWithClauses);
        }

        tok = tok(lexer, "'select', 'with' or table name expected");

        if (isWithKeyword(tok)) {
            parseWithClauses(lexer, model);
            tok = tok(lexer, "'select' or table name expected");
        }

        // [select]
        if (isSelectKeyword(tok)) {
            parseSelectClause(lexer, model);

            tok = optTok(lexer);

            if (tok != null && setOperations.contains(tok)) {
                tok = null;
            }

            if (tok == null) {
                QueryModel nestedModel = queryModelPool.next();
                nestedModel.setModelPosition(modelPosition);
                ExpressionNode func = expressionNodePool.next().of(ExpressionNode.FUNCTION, "long_sequence", 0, lexer.lastTokenPosition());
                func.paramCount = 1;
                func.rhs = expressionNodePool.next().of(ExpressionNode.CONSTANT, "1", 0, 0);
                nestedModel.setTableName(func);
                model.setSelectModelType(QueryModel.SELECT_MODEL_VIRTUAL);
                model.setNestedModel(nestedModel);
                lexer.unparse();
                return model;
            }
        } else {
            lexer.unparse();
            model.addBottomUpColumn(SqlUtil.nextColumn(queryColumnPool, expressionNodePool, "*", "*"));
        }

        QueryModel nestedModel = queryModelPool.next();
        nestedModel.setModelPosition(modelPosition);

        parseFromClause(lexer, nestedModel, model);
        if (nestedModel.getLimitHi() != null || nestedModel.getLimitLo() != null) {
            model.setLimit(nestedModel.getLimitLo(), nestedModel.getLimitHi());
            nestedModel.setLimit(null, null);
        }
        model.setSelectModelType(QueryModel.SELECT_MODEL_CHOOSE);
        model.setNestedModel(nestedModel);
        return model;
    }

    private void parseFromClause(GenericLexer lexer, QueryModel model, QueryModel masterModel) throws SqlException {
        CharSequence tok = expectTableNameOrSubQuery(lexer);
        // expect "(" in case of sub-query

        if (Chars.equals(tok, '(')) {
            model.setNestedModel(parseAsSubQueryAndExpectClosingBrace(lexer, masterModel.getWithClauses()));
            model.setNestedModelIsSubQuery(true);

            tok = optTok(lexer);

            // check if tok is not "where" - should be alias

            if (tok != null && tableAliasStop.excludes(tok)) {
                model.setAlias(literal(lexer, tok));
                tok = optTok(lexer);
            }

            // expect [timestamp(column)]

            ExpressionNode timestamp = parseTimestamp(lexer, tok);
            if (timestamp != null) {
                model.setTimestamp(timestamp);
                tok = optTok(lexer);
            }
        } else {

            lexer.unparse();
            parseSelectFrom(lexer, model, masterModel);

            tok = optTok(lexer);

            if (tok != null && tableAliasStop.excludes(tok)) {
                model.setAlias(literal(lexer, tok));
                tok = optTok(lexer);
            }

            // expect [timestamp(column)]

            ExpressionNode timestamp = parseTimestamp(lexer, tok);
            if (timestamp != null) {
                model.setTimestamp(timestamp);
                tok = optTok(lexer);
            }

            // expect [latest by]

            if (tok != null && isLatestKeyword(tok)) {
                parseLatestBy(lexer, model);
                tok = optTok(lexer);
            }
        }

        // expect multiple [[inner | outer | cross] join]

        int joinType;
        while (tok != null && (joinType = joinStartSet.get(tok)) != -1) {
            model.addJoinModel(parseJoin(lexer, tok, joinType, masterModel));
            tok = optTok(lexer);
        }

        // expect [where]

        if (tok != null && isWhereKeyword(tok)) {
            ExpressionNode expr = expr(lexer, model);
            if (expr != null) {
                model.setWhereClause(expr);
                tok = optTok(lexer);
            } else {
                throw SqlException.$((lexer.lastTokenPosition()), "empty where clause");
            }
        }

        // expect [sample by]

        if (tok != null && isSampleKeyword(tok)) {
            expectBy(lexer);
            model.setSampleBy(expectLiteral(lexer));
            tok = optTok(lexer);

            if (tok != null && isFillKeyword(tok)) {
                expectTok(lexer, '(');
                do {
                    final ExpressionNode fillNode = expr(lexer, model);
                    if (fillNode == null) {
                        throw SqlException.$(lexer.lastTokenPosition(), "'none', 'prev', 'mid', 'null' or number expected");
                    }
                    model.addSampleByFill(fillNode);
                    tok = tokIncludingLocalBrace(lexer, "',' or ')'");
                    if (Chars.equals(tok, ')')) {
                        break;
                    }
                    expectTok(tok, lexer.lastTokenPosition(), ',');
                } while (true);
                tok = optTok(lexer);
            }
        }

        //expect [group by]

        if (tok != null && isGroupKeyword(tok)) {
            expectBy(lexer);
            do {
                tokIncludingLocalBrace(lexer, "literal");
                lexer.unparse();
                ExpressionNode n = expr(lexer, model);
                if (n == null || (n.type != ExpressionNode.LITERAL && n.type != ExpressionNode.CONSTANT && n.type != ExpressionNode.FUNCTION && n.type != ExpressionNode.OPERATION)) {
                    throw SqlException.$(n == null ? lexer.lastTokenPosition() : n.position, "literal expected");
                }

                model.addGroupBy(n);

                tok = optTok(lexer);
            } while (tok != null && Chars.equals(tok, ','));
        }

        // expect [order by]

        if (tok != null && isOrderKeyword(tok)) {
            expectBy(lexer);
            do {
                tokIncludingLocalBrace(lexer, "literal");
                lexer.unparse();

                ExpressionNode n = expr(lexer, model);
                if (n == null || (n.type != ExpressionNode.LITERAL && n.type != ExpressionNode.CONSTANT)) {
                    throw SqlException.$(n == null ? lexer.lastTokenPosition() : n.position, "literal expected");
                }

                tok = optTok(lexer);

                if (tok != null && isDescKeyword(tok)) {

                    model.addOrderBy(n, QueryModel.ORDER_DIRECTION_DESCENDING);
                    tok = optTok(lexer);

                } else {

                    model.addOrderBy(n, QueryModel.ORDER_DIRECTION_ASCENDING);

                    if (tok != null && isAscKeyword(tok)) {
                        tok = optTok(lexer);
                    }
                }

                if (model.getOrderBy().size() >= MAX_ORDER_BY_COLUMNS) {
                    throw err(lexer, "Too many columns");
                }

            } while (tok != null && Chars.equals(tok, ','));
        }

        // expect [limit]
        if (tok != null && isLimitKeyword(tok)) {
            ExpressionNode lo = expr(lexer, model);
            ExpressionNode hi = null;

            tok = optTok(lexer);
            if (tok != null && Chars.equals(tok, ',')) {
                hi = expr(lexer, model);
            } else {
                lexer.unparse();
            }
            model.setLimit(lo, hi);
        } else {
            lexer.unparse();
        }
    }

    private ExecutionModel parseInsert(GenericLexer lexer) throws SqlException {
        expectTok(lexer, "into");

        final InsertModel model = insertModelPool.next();
        model.setTableName(expectLiteral(lexer));

        CharSequence tok = tok(lexer, "'(' or 'select'");

        if (Chars.equals(tok, '(')) {
            do {
                tok = tok(lexer, "column");
                if (Chars.equals(tok, ')')) {
                    throw err(lexer, "missing column name");
                }

                if (!model.addColumn(GenericLexer.immutableOf(GenericLexer.unquote(tok)), lexer.lastTokenPosition())) {
                    throw SqlException.position(lexer.lastTokenPosition()).put("duplicate column name: ").put(tok);
                }

            } while (Chars.equals((tok = tok(lexer, "','")), ','));

            expectTok(tok, lexer.lastTokenPosition(), ')');
            tok = optTok(lexer);
        }

        if (tok == null) {
            throw SqlException.$(lexer.getPosition(), "'select' or 'values' expected");
        }

        if (isSelectKeyword(tok)) {
            model.setSelectKeywordPosition(lexer.lastTokenPosition());
            lexer.unparse();
            final QueryModel queryModel = parseDml(lexer, null);
            model.setQueryModel(queryModel);
            return model;
        }

        if (isValuesKeyword(tok)) {
            expectTok(lexer, '(');

            do {
                ExpressionNode expr = expectExpr(lexer);
                if (Chars.equals(expr.token, ')')) {
                    throw err(lexer, "missing column value");
                }

                model.addColumnValue(expr);
            } while (Chars.equals((tok = tok(lexer, "','")), ','));

            expectTok(tok, lexer.lastTokenPosition(), ')');
            model.setEndOfValuesPosition(lexer.lastTokenPosition());

            return model;
        }
        throw err(lexer, "'select' or 'values' expected");
    }

    private ExecutionModel parseUpdate(GenericLexer lexer) throws SqlException {
        final UpdateModel model = updateModelObjectPool.next();
        QueryModel queryModel = queryModelPool.next();
        QueryModel queryNestedModel = queryModelPool.next();
        InsertModel insertModel = insertModelPool.next();
        queryModel.setNestedModel(queryNestedModel);
        model.setQueryModel(queryModel);
        model.setInsertModel(insertModel);
        queryModel.setSelectModelType(QueryModel.SELECT_MODEL_CHOOSE);

        model.setTableName(expectLiteral(lexer));

        expectTok(lexer, "set");
        CharSequence tok;

        do {
            tok = tok(lexer, "column");
            if (!model.addInsertColumn(GenericLexer.immutableOf(GenericLexer.unquote(tok)), lexer.lastTokenPosition())) {
                throw SqlException.position(lexer.lastTokenPosition()).put("duplicate column name: ").put(tok);
            }

            tok = tok(lexer, "=");
            if (tok != null && !Chars.equals(tok, '=')) {
                throw SqlException.position(lexer.lastTokenPosition()).put("expected '=' instead of: ").put(tok);
            }

            model.addInsertColumnValue(expectExpr(lexer));
        } while (Chars.equals((tok = tok(lexer, "','")), ','));

        if (tok != null && isWhereKeyword(tok)) {
            ExpressionNode expr = expr(lexer, queryNestedModel);
            if (expr != null) {
                queryNestedModel.setWhereClause(expr);
            } else {
                throw SqlException.$((lexer.lastTokenPosition()), "empty where clause");
            }
        } else {
            throw SqlException.$((lexer.lastTokenPosition()), "empty where clause");
        }

        return model;
    }

    private void parseLatestBy(GenericLexer lexer, QueryModel model) throws SqlException {
        expectBy(lexer);
        CharSequence tok;
        do {
            model.addLatestBy(expectLiteral(lexer));
            tok = SqlUtil.fetchNext(lexer);
        } while (Chars.equalsNc(tok, ','));

        if (tok != null) {
            lexer.unparse();
        }
    }

    private ExecutionModel parseRenameStatement(GenericLexer lexer) throws SqlException {
        expectTok(lexer, "table");
        RenameTableModel model = renameTableModelPool.next();
        ExpressionNode e = expectExpr(lexer);
        if (e.type != ExpressionNode.LITERAL && e.type != ExpressionNode.CONSTANT) {
            throw SqlException.$(e.position, "literal or constant expected");
        }
        model.setFrom(e);
        expectTok(lexer, "to");

        e = expectExpr(lexer);
        if (e.type != ExpressionNode.LITERAL && e.type != ExpressionNode.CONSTANT) {
            throw SqlException.$(e.position, "literal or constant expected");
        }
        model.setTo(e);
        return model;
    }

    private QueryModel parseJoin(GenericLexer lexer, CharSequence tok, int joinType, QueryModel parent) throws SqlException {
        QueryModel joinModel = queryModelPool.next();
        joinModel.setJoinType(joinType);
        joinModel.setJoinKeywordPosition(lexer.lastTokenPosition());

        if (!isJoinKeyword(tok) && !Chars.equals(tok, ',')) {
            expectTok(lexer, "join");
        }

        tok = expectTableNameOrSubQuery(lexer);

        if (Chars.equals(tok, '(')) {
            joinModel.setNestedModel(parseAsSubQueryAndExpectClosingBrace(lexer, parent.getWithClauses()));
        } else {
            lexer.unparse();
            parseSelectFrom(lexer, joinModel, parent);
        }

        tok = optTok(lexer);

        if (tok != null && tableAliasStop.excludes(tok)) {
            lexer.unparse();
            joinModel.setAlias(literal(lexer, optTok(lexer)));
        } else {
            lexer.unparse();
        }

        tok = optTok(lexer);

        if (joinType == QueryModel.JOIN_CROSS && tok != null && isOnKeyword(tok)) {
            throw SqlException.$(lexer.lastTokenPosition(), "Cross joins cannot have join clauses");
        }

        switch (joinType) {
            case QueryModel.JOIN_ASOF:
            case QueryModel.JOIN_LT:
            case QueryModel.JOIN_SPLICE:
                if (tok == null || !isOnKeyword(tok)) {
                    lexer.unparse();
                    break;
                }
                // intentional fall through
            case QueryModel.JOIN_INNER:
            case QueryModel.JOIN_OUTER:
                expectTok(lexer, tok, "on");
                try {
                    expressionParser.parseExpr(lexer, expressionTreeBuilder);
                    ExpressionNode expr;
                    switch (expressionTreeBuilder.size()) {
                        case 0:
                            throw SqlException.$(lexer.lastTokenPosition(), "Expression expected");
                        case 1:
                            expr = expressionTreeBuilder.poll();
                            if (expr.type == ExpressionNode.LITERAL) {
                                do {
                                    joinModel.addJoinColumn(expr);
                                } while ((expr = expressionTreeBuilder.poll()) != null);
                            } else {
                                joinModel.setJoinCriteria(rewriteCase(rewriteCount(expr)));
                            }
                            break;
                        default:
                            // this code handles "join on (a,b,c)", e.g. list of columns
                            while ((expr = expressionTreeBuilder.poll()) != null) {
                                if (expr.type != ExpressionNode.LITERAL) {
                                    throw SqlException.$(lexer.lastTokenPosition(), "Column name expected");
                                }
                                joinModel.addJoinColumn(expr);
                            }
                            break;
                    }
                } catch (SqlException e) {
                    expressionTreeBuilder.reset();
                    throw e;
                }
                break;
            default:
                lexer.unparse();
                break;
        }

        return joinModel;
    }

    private ExecutionModel parseSelect(GenericLexer lexer) throws SqlException {
        lexer.unparse();
        final QueryModel model = parseDml(lexer, null);
        final CharSequence tok = optTok(lexer);
        if (tok == null || Chars.equals(tok, ';')) {
            return model;
        }
        throw errUnexpected(lexer, tok);
    }

    private void parseSelectFrom(GenericLexer lexer, QueryModel model, QueryModel masterModel) throws SqlException {
        final ExpressionNode expr = expr(lexer, model);
        if (expr == null) {
            throw SqlException.position(lexer.lastTokenPosition()).put("table name expected");
        }
        CharSequence name = expr.token;

        switch (expr.type) {
            case ExpressionNode.LITERAL:
            case ExpressionNode.CONSTANT:
                final ExpressionNode literal = literal(name, expr.position);
                final WithClauseModel withClause = masterModel.getWithClause(name);
                if (withClause != null) {
                    model.setNestedModel(parseWith(lexer, withClause, masterModel.getWithClauses()));
                    model.setAlias(literal);
                } else {
                    model.setTableName(literal);
                }
                break;
            case ExpressionNode.FUNCTION:
                model.setTableName(expr);
                break;
            default:
                throw SqlException.$(expr.position, "function, literal or constant is expected");
        }
    }

    private int parseSymbolCapacity(GenericLexer lexer) throws SqlException {
        final int errorPosition = lexer.getPosition();
        final int symbolCapacity = expectInt(lexer);
        TableUtils.validateSymbolCapacity(errorPosition, symbolCapacity);
        return Numbers.ceilPow2(symbolCapacity);
    }

    private ExpressionNode parseTimestamp(GenericLexer lexer, CharSequence tok) throws SqlException {
        if (tok != null && isTimestampKeyword(tok)) {
            expectTok(lexer, '(');
            final ExpressionNode result = expectLiteral(lexer);
            tokIncludingLocalBrace(lexer, "')'");
            return result;
        }
        return null;
    }

    private QueryModel parseWith(GenericLexer lexer, WithClauseModel wcm, CharSequenceObjHashMap<WithClauseModel> withClauses) throws SqlException {
        QueryModel m = wcm.popModel();
        if (m != null) {
            return m;
        }

        final int pos = lexer.getPosition();
        final CharSequence unparsed = lexer.getUnparsed();
        lexer.goToPosition(wcm.getPosition(), null);
        // this will not throw exception because this is second pass over the same sub-query
        // we wouldn't be here is syntax was wrong
        m = parseAsSubQueryAndExpectClosingBrace(lexer, withClauses);
        lexer.goToPosition(pos, unparsed);
        return m;
    }

    private void parseWithClauses(GenericLexer lexer, QueryModel model) throws SqlException {
        do {
            ExpressionNode name = expectLiteral(lexer);

            if (model.getWithClause(name.token) != null) {
                throw SqlException.$(name.position, "duplicate name");
            }

            expectTok(lexer, "as");
            expectTok(lexer, '(');
            int lo = lexer.lastTokenPosition();
            WithClauseModel wcm = withClauseModelPool.next();
            wcm.of(lo + 1, parseAsSubQueryAndExpectClosingBrace(lexer, model.getWithClauses()));
            model.addWithClause(name.token, wcm);

            CharSequence tok = optTok(lexer);
            if (tok == null || !Chars.equals(tok, ',')) {
                lexer.unparse();
                break;
            }
        } while (true);
    }

    private ExpressionNode rewriteCase(ExpressionNode parent) throws SqlException {
        traversalAlgo.traverse(parent, rewriteCase0Ref);
        return parent;
    }

    private void rewriteCase0(ExpressionNode node) {
        if (node.type == ExpressionNode.FUNCTION && isCaseKeyword(node.token)) {
            tempExprNodes.clear();
            ExpressionNode literal = null;
            ExpressionNode elseExpr;
            boolean convertToSwitch = true;
            final int paramCount = node.paramCount;

            if (node.paramCount == 2) {
                // special case, typically something like
                // case value else expression end
                // this can be simplified to "expression" only

                ExpressionNode that = node.rhs;
                node.of(that.type, that.token, that.precedence, that.position);
                node.paramCount = that.paramCount;
                if (that.paramCount == 2) {
                    node.lhs = that.lhs;
                    node.rhs = that.rhs;
                } else {
                    node.args.clear();
                    node.args.addAll(that.args);
                }
                return;
            }
            final int lim;
            if ((paramCount & 1) == 0) {
                elseExpr = node.args.getQuick(0);
                lim = 0;
            } else {
                elseExpr = null;
                lim = -1;
            }

            // agrs are in inverted order, hence last list item is the first arg
            ExpressionNode first = node.args.getQuick(paramCount - 1);
            if (first.token != null) {
                // simple case of 'case' :) e.g.
                // case x
                //   when 1 then 'A'
                //   ...
                node.token = "switch";
                return;
            }

            for (int i = paramCount - 2; i > lim; i--) {
                if ((i & 1) == 1) {
                    // this is "then" clause, copy it as as
                    tempExprNodes.add(node.args.getQuick(i));
                    continue;
                }
                ExpressionNode where = node.args.getQuick(i);
                if (where.type == ExpressionNode.OPERATION && where.token.charAt(0) == '=') {
                    ExpressionNode thisConstant;
                    ExpressionNode thisLiteral;
                    if (where.lhs.type == ExpressionNode.CONSTANT && where.rhs.type == ExpressionNode.LITERAL) {
                        thisConstant = where.lhs;
                        thisLiteral = where.rhs;
                    } else if (where.lhs.type == ExpressionNode.LITERAL && where.rhs.type == ExpressionNode.CONSTANT) {
                        thisConstant = where.rhs;
                        thisLiteral = where.lhs;
                    } else {
                        convertToSwitch = false;
                        // not supported
                        break;
                    }

                    if (literal == null) {
                        literal = thisLiteral;
                        tempExprNodes.add(thisConstant);
                    } else if (Chars.equals(literal.token, thisLiteral.token)) {
                        tempExprNodes.add(thisConstant);
                    } else {
                        convertToSwitch = false;
                        // not supported
                        break;
                    }
                } else {
                    convertToSwitch = false;
                    // not supported
                    break;
                }
            }

            if (convertToSwitch) {
                int n = tempExprNodes.size();
                node.token = "switch";
                node.args.clear();
                node.args.add(elseExpr);
                for (int i = n - 1; i > -1; i--) {
                    node.args.add(tempExprNodes.getQuick(i));
                }
                node.args.add(literal);
                node.paramCount = n + 2;
            } else {
                node.args.remove(paramCount - 1);
                node.paramCount = paramCount - 1;
            }
        }
    }

    private ExpressionNode rewriteCount(ExpressionNode parent) throws SqlException {
        traversalAlgo.traverse(parent, rewriteCount0Ref);
        return parent;
    }

    /**
     * Rewrites count(*) expressions to count().
     *
     * @param node expression node, provided by tree walking algo
     */
    private void rewriteCount0(ExpressionNode node) {
        if (node.type == ExpressionNode.FUNCTION && isCountKeyword(node.token)) {
            if (node.paramCount == 1) {
                // special case, typically something like
                // case value else expression end
                // this can be simplified to "expression" only

                ExpressionNode that = node.rhs;
                if (Chars.equals(that.token, '*')) {
                    if (that.rhs == null && node.lhs == null) {
                        that.paramCount = 0;
                        node.rhs = null;
                        node.paramCount = 0;
                    }
                }
            }
        }
    }

    private int toColumnType(GenericLexer lexer, CharSequence tok) throws SqlException {
        final int type = ColumnType.columnTypeOf(tok);
        if (type == -1) {
            throw SqlException.$(lexer.lastTokenPosition(), "unsupported column type: ").put(tok);
        }
        return type;
    }

    private @NotNull CharSequence tok(GenericLexer lexer, String expectedList) throws SqlException {
        final int pos = lexer.getPosition();
        CharSequence tok = optTok(lexer);
        if (tok == null) {
            throw SqlException.position(pos).put(expectedList).put(" expected");
        }
        return tok;
    }

    private @NotNull CharSequence tokIncludingLocalBrace(GenericLexer lexer, String expectedList) throws SqlException {
        final int pos = lexer.getPosition();
        final CharSequence tok = SqlUtil.fetchNext(lexer);
        if (tok == null) {
            throw SqlException.position(pos).put(expectedList).put(" expected");
        }
        return tok;
    }

    private void validateLiteral(int pos, CharSequence tok) throws SqlException {
        switch (tok.charAt(0)) {
            case '(':
            case ')':
            case ',':
            case '`':
//            case '"':
            case '\'':
                throw SqlException.position(pos).put("literal expected");
            default:
                break;

        }
    }
}
