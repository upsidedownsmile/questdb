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

package io.questdb.griffin.engine.functions.lt;

import io.questdb.cairo.CairoConfiguration;
import io.questdb.cairo.CairoException;
import io.questdb.cairo.TableUtils;
import io.questdb.cairo.sql.Function;
import io.questdb.cairo.sql.Record;
import io.questdb.griffin.AbstractBooleanFunctionFactory;
import io.questdb.griffin.FunctionFactory;
import io.questdb.griffin.SqlException;
import io.questdb.griffin.engine.functions.BinaryFunction;
import io.questdb.griffin.engine.functions.BooleanFunction;
import io.questdb.griffin.engine.functions.UnaryFunction;
import io.questdb.std.NumericException;
import io.questdb.std.ObjList;
import io.questdb.std.microtime.TimestampFormatUtils;

public class LtTimestampStringFunctionFactory extends AbstractBooleanFunctionFactory implements FunctionFactory {
    @Override
    public String getSignature() {
        return "<(NS)";
    }

    @Override
    public Function newInstance(ObjList<Function> args, int position, CairoConfiguration configuration) throws SqlException {
        Function strFn = args.getQuick(1);
        if (strFn.isConstant()) {
            try {
                long value = TableUtils.toTimestampOrException(strFn.getStr(null));
                return new LtTimestampVVFunction(position, args.getQuick(0), value, isNegated);
            } catch (CairoException e) {
                throw SqlException.$(position, e.getFlyweightMessage());
            }

        }
        return new LtTimestampStringFunction(position, args.getQuick(0), args.getQuick(1), isNegated);
    }

    private static class LtTimestampStringFunction extends BooleanFunction implements BinaryFunction {
        private final boolean isNegated;
        private final Function left;
        private final Function right;

        public LtTimestampStringFunction(int position, Function left, Function right, boolean isNegated) {
            super(position);
            this.left = left;
            this.right = right;
            this.isNegated = isNegated;
        }

        @Override
        public boolean getBool(Record rec) {
            long value;
            final CharSequence str = right.getStr(rec);
            try {
                value = TimestampFormatUtils.parseTimestamp(str);
            } catch (NumericException e) {
                try {
                    value = TimestampFormatUtils.parseDateTime(str);
                } catch (NumericException numericException) {
                    return false;
                }
            }
            return isNegated == (left.getTimestamp(rec) >= value);
        }

        @Override
        public Function getLeft() {
            return left;
        }

        @Override
        public Function getRight() {
            return right;
        }
    }

    private static class LtTimestampVVFunction extends BooleanFunction implements UnaryFunction {
        private final boolean isNegated;
        private final Function left;
        private final long right;

        public LtTimestampVVFunction(int position, Function left, long right, boolean isNegated) {
            super(position);
            this.left = left;
            this.right = right;
            this.isNegated = isNegated;
        }

        @Override
        public Function getArg() {
            return left;
        }

        @Override
        public boolean getBool(Record rec) {
            return isNegated == (left.getTimestamp(rec) >= right);
        }
    }
}
