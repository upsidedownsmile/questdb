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

package io.questdb.griffin.model;

import io.questdb.std.CharSequenceHashSet;
import io.questdb.std.IntList;
import io.questdb.std.Mutable;
import io.questdb.std.ObjList;
import io.questdb.std.ObjectFactory;
import io.questdb.std.Sinkable;
import io.questdb.std.str.CharSink;

public class UpdateModel implements Mutable, ExecutionModel, Sinkable {
    public static final ObjectFactory<UpdateModel> FACTORY = UpdateModel::new;

    private final CharSequenceHashSet columnSet = new CharSequenceHashSet();

    private final ObjList<ExpressionNode> columnValues = new ObjList<>();

    private final IntList columnPositions = new IntList();

    private ExpressionNode tableName;

    private QueryModel queryModel;

    private UpdateModel() {
    }

    public void setTableName(ExpressionNode tableName) {
        this.tableName = tableName;
    }

    public boolean addColumn(CharSequence columnName, int columnPosition) {
        if (columnSet.add(columnName)) {
            columnPositions.add(columnPosition);
            return true;
        }
        return false;
    }

    public QueryModel getQueryModel() {
        return queryModel;
    }

    public void setQueryModel(QueryModel queryModel) {
        this.queryModel = queryModel;
    }


    public void addColumnValue(ExpressionNode value) {
        columnValues.add(value);
    }

    @Override
    public void clear() {
        this.tableName = null;
        this.queryModel = null;
        this.columnSet.clear();
        this.columnPositions.clear();
        this.columnValues.clear();
    }

    @Override
    public int getModelType() {
        return UPDATE;
    }

    @Override
    public void toSink(final CharSink sink) {
        sink.put("update ").put(tableName.token).put(" set ");
        int n = columnSet.size();
        for(int idx = 0; idx < n; idx++) {
            sink.put(columnSet.get(idx))
                .put(" = ")
                .put(columnValues.get(idx));

            if(idx != n-1) {
                sink.put(", ");
            }
        }
        sink.put(" where ");
        queryModel.getWhereClause().toSink(sink);
    }
}
