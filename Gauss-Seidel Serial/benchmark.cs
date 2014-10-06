using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;

namespace Hill_Cipher.Research
{
    class benchmark
    {
        public benchmark()
        {
            stopWatch = new Stopwatch();
        }

        Stopwatch stopWatch;

        public void startBenchmark()
        {
            stopWatch.Reset();
            stopWatch.Start();
        }

        public string getElapsedTime()
        {
            TimeSpan ts = stopWatch.Elapsed;
            return format(ts);
        }

        public string getBenchmarkResult()
        {
            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            return format(ts);
        }

        private string format(TimeSpan ts)
        {
            return (String.Format("{0:00}:{1:00}:{2:00}.{3:000}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds));
        }
    }
}
