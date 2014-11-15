using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;

namespace Gauss_Seidel_Sequential
{
    class benchmark
    {
        public benchmark()
        {
            stopWatch = new Stopwatch();
        }

        Stopwatch stopWatch;

        public void start()
        {
            stopWatch.Reset();
            stopWatch.Start();
        }

        public void pause()
        {
            stopWatch.Stop();
        }

        public void resume()
        {
            stopWatch.Start();
        }

        public string getElapsedTime()
        {
            TimeSpan ts = stopWatch.Elapsed;
            return format(ts);
        }

        public string getResult()
        {
            stopWatch.Stop();
            TimeSpan ts = stopWatch.Elapsed;
            return format(ts);
        }

        private string format(TimeSpan ts)
        {
            return (String.Format("{0:00}:{1:00}:{2:00}.{3:000}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds));
        }

        public double getElapsedSeconds()
        {
            TimeSpan ts = stopWatch.Elapsed;
            return ts.TotalSeconds;
        }
    }
}
