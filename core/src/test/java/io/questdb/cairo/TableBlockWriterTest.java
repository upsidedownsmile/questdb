package io.questdb.cairo;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.LockSupport;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import io.questdb.cairo.TableBlockWriter.TableBlockWriterJob;
import io.questdb.cairo.security.AllowAllCairoSecurityContext;
import io.questdb.cairo.sql.PageFrame;
import io.questdb.cairo.sql.RecordCursor;
import io.questdb.cairo.sql.RecordCursorFactory;
import io.questdb.griffin.AbstractGriffinTest;
import io.questdb.griffin.CompiledQuery;
import io.questdb.griffin.SqlException;
import io.questdb.log.Log;
import io.questdb.log.LogFactory;
import io.questdb.mp.SOCountDownLatch;
import io.questdb.test.tools.TestUtils;
import io.questdb.test.tools.TestUtils.LeakProneCode;

public class TableBlockWriterTest extends AbstractGriffinTest {
    private static final Log LOG = LogFactory.getLog(TableBlockWriterTest.class);

    @Test
    public void testSimple() throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_long(-55, 9009, 2) l FROM long_sequence(500)" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, l LONG) TIMESTAMP(ts);", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testSimpleResumeBlock() throws Exception {
        runTest(false, () -> {
            int nConsecutiveRows = 50;
            long tsStart = 0;
            long tsInc = 1000000000;
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT" +
                    " rnd_long(100,200,2) j," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts" +
                    " from long_sequence(" + nConsecutiveRows + ")" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile(
                    "CREATE TABLE dest (j LONG, ts TIMESTAMP) TIMESTAMP(ts);",
                    sqlExecutionContext);
            replicateTable("source", "dest");
            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            tsStart += nConsecutiveRows * tsInc;
            compiler.compile("INSERT INTO source(j, ts) " +
                    "SELECT" +
                    " rnd_long(100,200,2) j," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts" +
                    " from long_sequence(" + nConsecutiveRows + ")" +
                    ";",
                    sqlExecutionContext);
            expected = select("SELECT * FROM source");
            replicateTable("source", "dest", nConsecutiveRows);
            actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testSimpleResumeBlockWithRetry() throws Exception {
        testSimpleResumeBlockWithRetry(true);
    }

    @Test
    public void testSimpleResumeBlockWithRetryAndNoCancel() throws Exception {
        testSimpleResumeBlockWithRetry(false);
    }

    private void testSimpleResumeBlockWithRetry(boolean cancel) throws Exception {
        runTest(false, () -> {
            int nConsecutiveRows = 10;
            long tsStart = 0;
            long tsInc = 1000000000;
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT" +
                    " rnd_long(100,200,2) j," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts" +
                    " from long_sequence(" + nConsecutiveRows + ")" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile(
                    "CREATE TABLE dest (j LONG, ts TIMESTAMP) TIMESTAMP(ts);",
                    sqlExecutionContext);
            replicateTable("source", "dest");
            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            tsStart += nConsecutiveRows * tsInc;
            compiler.compile("INSERT INTO source(j, ts) " +
                    "SELECT" +
                    " rnd_long(100,200,2) j," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts" +
                    " from long_sequence(" + nConsecutiveRows + ")" +
                    ";",
                    sqlExecutionContext);
            replicateTable("source", "dest", nConsecutiveRows, false, Long.MAX_VALUE, cancel);
            actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            replicateTable("source", "dest", nConsecutiveRows, true, Long.MAX_VALUE, false);
            actual = select("SELECT * FROM dest");
            expected = select("SELECT * FROM source");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testPartitioned() throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_long(-55, 9009, 2) l FROM long_sequence(500)" +
                    ") TIMESTAMP (ts) PARTITION BY DAY;",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, l LONG) TIMESTAMP(ts) PARTITION BY DAY;", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testNoTimestamp() throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT rnd_long(-55, 9009, 2) l FROM long_sequence(500)" +
                    ");",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (l LONG);", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testString1() throws Exception {
        testString(false);
    }

    @Test
    public void testString2() throws Exception {
        testString(true);
    }

    private void testString(boolean endsWithNull) throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_str(5,10,2) s FROM long_sequence(300)" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            if (endsWithNull) {
                compiler.compile("INSERT INTO source (ts, s) SELECT ts+1, null FROM (" +
                        "SELECT ts, s FROM source ORDER BY ts DESC LIMIT 1" +
                        ")",
                        sqlExecutionContext);
            } else {
                compiler.compile("INSERT INTO source (ts, s) SELECT ts+1, 'ABC' FROM (" +
                        "SELECT ts, s FROM source ORDER BY ts DESC LIMIT 1" +
                        ")",
                        sqlExecutionContext);
            }
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, s STRING) TIMESTAMP(ts);", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testBinary1() throws Exception {
        testBinary(false);
    }

    @Test
    public void testBinary2() throws Exception {
        testBinary(true);
    }

    private void testBinary(boolean endsWithNull) throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_bin(10, 20, 2) bin FROM long_sequence(500)" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            if (endsWithNull) {
                compiler.compile("INSERT INTO source (ts) SELECT ts+1 FROM (" +
                        "SELECT ts, bin FROM source ORDER BY ts DESC LIMIT 1" +
                        ")",
                        sqlExecutionContext);
            } else {
                compiler.compile("INSERT INTO source (ts, bin) SELECT ts+1, rnd_bin() FROM (" +
                        "SELECT ts, bin FROM source ORDER BY ts DESC LIMIT 1" +
                        ")",
                        sqlExecutionContext);
            }
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, bin BINARY) TIMESTAMP(ts);", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testSymbol() throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_symbol(60,2,16,2) sym FROM long_sequence(500)" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, sym SYMBOL) TIMESTAMP(ts);", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAllTypes() throws Exception {
        runTest(true, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(0, 1000000000) ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(1000)" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile(
                    "CREATE TABLE dest (ch CHAR, ll LONG256, a1 INT, a INT, b BOOLEAN, c STRING, d DOUBLE, e FLOAT, f SHORT, f1 SHORT, g DATE, h TIMESTAMP, i SYMBOL, j LONG, j1 LONG, ts TIMESTAMP, l BYTE, m BINARY) TIMESTAMP(ts);",
                    sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAllTypesPartitioned() throws Exception {
        runTest(true, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(0, 1000000000) ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(1000)" +
                    ") TIMESTAMP (ts) PARTITION BY DAY;",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile(
                    "CREATE TABLE dest (ch CHAR, ll LONG256, a1 INT, a INT, b BOOLEAN, c STRING, d DOUBLE, e FLOAT, f SHORT, f1 SHORT, g DATE, h TIMESTAMP, i SYMBOL, j LONG, j1 LONG, ts TIMESTAMP, l BYTE, m BINARY) TIMESTAMP(ts) PARTITION BY DAY;",
                    sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAllTypesResumeBlock() throws Exception {
        testAllTypesResumeBlock(Long.MAX_VALUE);
    }

    @Test
    public void testAllTypesResumeBlockFragmentedFrames() throws Exception {
        testAllTypesResumeBlock(4);
    }

    public void testAllTypesResumeBlock(long maxRowsPerFrame) throws Exception {
        runTest(true, () -> {
            int nConsecutiveRows = 50;
            long tsStart = 0;
            long tsInc = 1000000000;
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(" + nConsecutiveRows + ")" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile(
                    "CREATE TABLE dest (ch CHAR, ll LONG256, a1 INT, a INT, b BOOLEAN, c STRING, d DOUBLE, e FLOAT, f SHORT, f1 SHORT, g DATE, h TIMESTAMP, i SYMBOL, j LONG, j1 LONG, ts TIMESTAMP, l BYTE, m BINARY) TIMESTAMP(ts);",
                    sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            tsStart += nConsecutiveRows * tsInc;
            compiler.compile("INSERT INTO source(ch, ll, a1, a, b, c, d, e, f, f1, g, h, i, j, j1, ts, l, m) " +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(" + nConsecutiveRows + ")" +
                    ";",
                    sqlExecutionContext);
            expected = select("SELECT * FROM source");
            replicateTable("source", "dest", nConsecutiveRows);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAllTypesPartitionedResumeBlock() throws Exception {
        testAllTypesPartitionedResumeBlock(Long.MAX_VALUE);
    }

    @Test
    public void testAllTypesPartitionedResumeBlockFragmentedFrames() throws Exception {
        testAllTypesPartitionedResumeBlock(3);
    }

    private void testAllTypesPartitionedResumeBlock(long maxRowsPerFrame) throws Exception {
        runTest(true, () -> {
            long nConsecuriveRows = 50;
            long tsStart = 0;
            long tsInc = 1000000000;
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(" + nConsecuriveRows + ")" +
                    ") TIMESTAMP (ts) PARTITION BY DAY;",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile(
                    "CREATE TABLE dest (ch CHAR, ll LONG256, a1 INT, a INT, b BOOLEAN, c STRING, d DOUBLE, e FLOAT, f SHORT, f1 SHORT, g DATE, h TIMESTAMP, i SYMBOL, j LONG, j1 LONG, ts TIMESTAMP, l BYTE, m BINARY) TIMESTAMP(ts) PARTITION BY DAY;",
                    sqlExecutionContext);
            replicateTable("source", "dest", 0, true, maxRowsPerFrame, false);

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            tsStart += nConsecuriveRows * tsInc;
            compiler.compile("INSERT INTO source(ch, ll, a1, a, b, c, d, e, f, f1, g, h, i, j, j1, ts, l, m) " +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(" + (nConsecuriveRows + 2) + ")" +
                    ";",
                    sqlExecutionContext);
            expected = select("SELECT * FROM source");
            replicateTable("source", "dest", nConsecuriveRows, true, maxRowsPerFrame, false);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAllTypesResumeBlockWithRetry() throws Exception {
        testAllTypesResumeBlockWithRetry(true);
    }

    @Test
    public void testAllTypesResumeBlockWithRetryAndNoCancel() throws Exception {
        testAllTypesResumeBlockWithRetry(false);
    }

    private void testAllTypesResumeBlockWithRetry(boolean cancel) throws Exception {
        runTest(true, () -> {
            int nConsecuriveRows = 50;
            long tsStart = 0;
            long tsInc = 1000000000;
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(" + nConsecuriveRows + ")" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile(
                    "CREATE TABLE dest (ch CHAR, ll LONG256, a1 INT, a INT, b BOOLEAN, c STRING, d DOUBLE, e FLOAT, f SHORT, f1 SHORT, g DATE, h TIMESTAMP, i SYMBOL, j LONG, j1 LONG, ts TIMESTAMP, l BYTE, m BINARY) TIMESTAMP(ts);",
                    sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            tsStart += nConsecuriveRows * tsInc;
            compiler.compile("INSERT INTO source(ch, ll, a1, a, b, c, d, e, f, f1, g, h, i, j, j1, ts, l, m) " +
                    "SELECT" +
                    " rnd_char() ch," +
                    " rnd_long256() ll," +
                    " rnd_int() a1," +
                    " rnd_int(0, 30, 2) a," +
                    " rnd_boolean() b," +
                    " rnd_str(3,3,2) c," +
                    " rnd_double(2) d," +
                    " rnd_float(2) e," +
                    " rnd_short(10,1024) f," +
                    " rnd_short() f1," +
                    " rnd_date(to_date('2015', 'yyyy'), to_date('2016', 'yyyy'), 2) g," +
                    " rnd_timestamp(to_timestamp('2015', 'yyyy'), to_timestamp('2016', 'yyyy'), 2) h," +
                    " rnd_symbol(4,4,4,2) i," +
                    " rnd_long(100,200,2) j," +
                    " rnd_long() j1," +
                    " timestamp_sequence(" + tsStart + ", " + tsInc + ") ts," +
                    " rnd_byte(2,50) l," +
                    " rnd_bin(10, 20, 2) m" +
                    " from long_sequence(" + nConsecuriveRows + ")" +
                    ";",
                    sqlExecutionContext);

            replicateTable("source", "dest", nConsecuriveRows, false, Long.MAX_VALUE, cancel);
            actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            expected = select("SELECT * FROM source");
            replicateTable("source", "dest", nConsecuriveRows, true, Long.MAX_VALUE, false);
            actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAddColumn1() throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_long(-55, 9009, 2) l FROM long_sequence(5)" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            compiler.compile("ALTER TABLE source ADD COLUMN str STRING",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, l LONG, str STRING) TIMESTAMP(ts);", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAddColumn2() throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_long(-55, 9009, 2) l FROM long_sequence(5)" +
                    ") TIMESTAMP (ts);",
                    sqlExecutionContext);
            compiler.compile("ALTER TABLE source ADD COLUMN str STRING",
                    sqlExecutionContext);
            compiler.compile("INSERT INTO source(ts, l, str) " +
                    "SELECT" +
                    " timestamp_sequence(5000000000, 500000000) ts," +
                    " rnd_long(-55, 9009, 2) l," +
                    " rnd_str(3,3,2) str" +
                    " from long_sequence(5)" +
                    ";",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, l LONG, str STRING) TIMESTAMP(ts);", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    @Test
    public void testAddColumnPartitioned() throws Exception {
        runTest(false, () -> {
            compiler.compile("CREATE TABLE source AS (" +
                    "SELECT timestamp_sequence(0, 1000000000) ts, rnd_long(-55, 9009, 2) l FROM long_sequence(200)" +
                    ") TIMESTAMP (ts) PARTITION BY DAY;",
                    sqlExecutionContext);
            compiler.compile("ALTER TABLE source ADD COLUMN str STRING",
                    sqlExecutionContext);
            compiler.compile("INSERT INTO source(ts, l, str) " +
                    "SELECT" +
                    " timestamp_sequence(400000000000, 500000000) ts," +
                    " rnd_long(-55, 9009, 2) l," +
                    " rnd_str(3,3,2) str" +
                    " from long_sequence(250)" +
                    ";",
                    sqlExecutionContext);
            String expected = select("SELECT * FROM source");

            compiler.compile("CREATE TABLE dest (ts TIMESTAMP, l LONG, str STRING) TIMESTAMP(ts) PARTITION BY DAY;", sqlExecutionContext);
            replicateTable("source", "dest");

            String actual = select("SELECT * FROM dest");
            Assert.assertEquals(expected, actual);

            engine.releaseInactive();
        });
    }

    private final AtomicInteger busyCount = new AtomicInteger();

    private void runTest(boolean threaded, LeakProneCode runnable) throws Exception {
        final SOCountDownLatch running = new SOCountDownLatch(1);
        final SOCountDownLatch finished = new SOCountDownLatch(1);
        Thread t = null;
        if (threaded) {
            busyCount.set(0);
            final TableBlockWriterJob job = new TableBlockWriterJob(engine.getMessageBus());
            t = new Thread() {
                @Override
                public void run() {
                    while (running.getCount() > 0) {
                        boolean busy = job.run(0);
                        if (busy) {
                            busyCount.incrementAndGet();
                        }
                        Thread.yield();
                    }
                    finished.countDown();
                }
            };
            t.start();
        } else {
            busyCount.set(1);
        }
        TestUtils.assertMemoryLeak(runnable);
        if (null != t) {
            running.countDown();
            finished.await();
        }
    }

    private void replicateTable(String sourceTableName, String destTableName) {
        replicateTable(sourceTableName, destTableName, 0);
    }

    private void replicateTable(String sourceTableName, String destTableName, long nFirstRow) {
        replicateTable(sourceTableName, destTableName, nFirstRow, true, Long.MAX_VALUE, false);
    }

    private void replicateTable(String sourceTableName, String destTableName, long nFirstRow, boolean commit, long maxRowsPerFrame, boolean cancel) {
        LOG.info().$("Replicating table from row ").$(nFirstRow).$();
        try (TableReplicationRecordCursorFactory factory = createReplicatingRecordCursorFactory(sourceTableName, maxRowsPerFrame);
                TableReader reader = engine.getReader(AllowAllCairoSecurityContext.INSTANCE, sourceTableName);
                TableWriter writer = engine.getWriter(sqlExecutionContext.getCairoSecurityContext(), destTableName);) {

            final int columnCount = writer.getMetadata().getColumnCount();
            int nFrames = 0;
            int timestampColumnIndex = reader.getMetadata().getTimestampIndex();
            TablePageFrameCursor cursor = factory.getPageFrameCursorFrom(sqlExecutionContext, timestampColumnIndex, nFirstRow);
            PageFrame frame;
            while ((frame = cursor.next()) != null) {
                TableBlockWriter blockWriter = writer.newBlock();
                long firstTimestamp = frame.getFirstTimestamp();
                long blockNRows = frame.getPageValueCount(0);
                LOG.info().$("Replicating frame from ").$ts(firstTimestamp).$();
                blockWriter.startPageFrame(firstTimestamp, blockNRows);

                for (int columnIndex = 0; columnIndex < columnCount; columnIndex++) {
                    int columnType = writer.getMetadata().getColumnType(columnIndex);
                    if (columnType == ColumnType.SYMBOL) {
                        SymbolMapReader symReader = reader.getSymbolMapReader(columnIndex);
                        int nSourceSymbols = symReader.size();
                        int nDestinationSymbols = writer.getSymbolMapWriter(columnIndex).getSymbolCount();

                        if (nSourceSymbols > nDestinationSymbols) {
                            long address = symReader.symbolCharsAddressOf(nDestinationSymbols);
                            long addressHi = symReader.symbolCharsAddressOf(nSourceSymbols);
                            blockWriter.appendSymbolCharsBlock(columnIndex, addressHi - address, address);
                        }
                    }
                }

                for (int columnIndex = 0; columnIndex < columnCount; columnIndex++) {
                    long pageAddress = frame.getPageAddress(columnIndex);
                    long blockLength = frame.getPageLength(columnIndex);
                    blockWriter.appendPageFrameColumn(columnIndex, blockLength, pageAddress);
                }
                nFrames++;
                if (commit) {
                    while (busyCount.get() == 0) {
                        LockSupport.parkNanos(0);
                    }
                    assert !cancel;
                    blockWriter.commit();
                } else {
                    if (cancel) {
                        blockWriter.cancel();
                    }
                    break;
                }
            }
            LOG.info().$("Replication finished in ").$(nFrames).$(" frames per row").$();
        }
    }

    private String select(CharSequence selectSql) throws SqlException {
        sink.clear();
        CompiledQuery query = compiler.compile(selectSql, sqlExecutionContext);
        try (RecordCursorFactory factory = query.getRecordCursorFactory(); RecordCursor cursor = factory.getCursor(sqlExecutionContext);) {
            printer.print(cursor, factory.getMetadata(), true);
        }
        return sink.toString();
    }

    private TableReplicationRecordCursorFactory createReplicatingRecordCursorFactory(String tableName, long maxRowsPerFrame) {
        return new TableReplicationRecordCursorFactory(engine, tableName, maxRowsPerFrame);
    }

    @BeforeClass
    public static void setUp3() {
        sqlExecutionContext.getRandom().reset(0, 1);
    }
}
