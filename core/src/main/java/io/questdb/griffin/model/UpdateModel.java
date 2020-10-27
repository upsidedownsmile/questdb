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
import io.questdb.std.Mutable;
import io.questdb.std.ObjectFactory;
import io.questdb.std.Sinkable;
import io.questdb.std.str.CharSink;

public class UpdateModel implements ExecutionModel, Mutable, Sinkable {
    public static final ObjectFactory<UpdateModel> FACTORY = UpdateModel::new;

    private QueryModel queryModel;

    private InsertModel insertModel;

    private UpdateModel() {
    }

    @Override
    public void clear() {
        this.queryModel = null;
        this.insertModel = null;
    }

    public QueryModel getQueryModel() {
        return queryModel;
    }

    public void setQueryModel(final QueryModel queryModel) {
        this.queryModel = queryModel;
    }

    public InsertModel getInsertModel() {
        return insertModel;
    }

    public void setInsertModel(final InsertModel insertModel) {
        this.insertModel = insertModel;
    }

    public void setTableName(ExpressionNode tableName) {
        queryModel.getNestedModel().setTableName(tableName);
        insertModel.setTableName(tableName);
    }

    public ExpressionNode getTableName() {
        return insertModel.getTableName();
    }

    public boolean addInsertColumn(CharSequence columnName, int columnPosition) {
        return insertModel.addColumn(columnName, columnPosition);
    }

    public void addInsertColumnValue(ExpressionNode value) {
        insertModel.addColumnValue(value);
    }

    public void setWhereClause(ExpressionNode whereClause) {
        queryModel.setWhereClause(whereClause);
    }

    public CharSequenceHashSet getUpdatedColumns() {
        return insertModel.getColumnSet();
    }

    @Override
    public void toSink(final CharSink sink) {

    }

    @Override
    public int getModelType() {
        return UPDATE;
    }
}
