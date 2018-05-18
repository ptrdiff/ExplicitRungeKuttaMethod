using System;
using MathNet.Numerics.LinearAlgebra;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ExplicitRungeKuttaMethod
{


    class ChartForm : Form
    {
        Chart myChart = new Chart();
        public ChartForm()
        {
            myChart.Parent = this;
            myChart.Dock = DockStyle.Fill;
            First();
            Second();
            Third();
        }
        private void First()
        {
            myChart.ChartAreas.Add(new ChartArea("First"));
            AddFullErrOfK("First");
            AddLine("First");
        }

        private void Second()
        {
            myChart.ChartAreas.Add(new ChartArea("Second"));
            AddFullErrOfInd("Second");
        }

        private void Third()
        {
            myChart.ChartAreas.Add(new ChartArea("Third"));
        }

        private void AddFullErrOfK(string chartArea)
        {
            Series mySeriesOfPoint = new Series("FullErrOfK")
            {
                ChartType = SeriesChartType.Line,
                ChartArea = chartArea
            };
            for (int k = 5; k <= 10; ++k)
            {
                var x = Math.Pow(2.0, -k);
                var y = (ERKMethod.RealFunction(5) - ERKMethod.MyMethod(x)).L2Norm();
                mySeriesOfPoint.Points.AddXY(y, x);
            }
            myChart.Series.Add(mySeriesOfPoint);
        }
        private void AddLine(string chartArea)
        {
            Series mySeriesOfPoint = new Series("Line")
            {
                ChartType = SeriesChartType.Line,
                ChartArea = chartArea
            };
            for (int k = 5; k <= 10; ++k)
            {
                var x = Math.Pow(2.0, -k);
                var y = 2 * x;
                mySeriesOfPoint.Points.AddXY(x, y);
            }
            myChart.Series.Add(mySeriesOfPoint);
        }
        private void AddFullErrOfInd(string chartArea)
        {
            Series mySeriesOfPoint = new Series("FullErrOfInd")
            {
                ChartType = SeriesChartType.Line,
                ChartArea = chartArea
            };

            for (int i = 1; i <= 5; ++i)
            {
                double h = ERKMethod.RungeRuleGlobal(0.00001, b: i);
                var y = (ERKMethod.RealFunction(i) - ERKMethod.MyMethod(h, b: i)).L2Norm();
                mySeriesOfPoint.Points.AddXY(i, y);
            }
            myChart.Series.Add(mySeriesOfPoint);
        }

    }

    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine(ERKMethod.RealFunction(5));
            Console.WriteLine(ERKMethod.MyMethod(ERKMethod.RungeRuleGlobal(0.00001)));
            ChartForm window = new ChartForm();
            window.ShowDialog();
        }
    }

    class ERKMethod
    {
        private static double A = 3.0;
        private static double B = 1.5;
        private static double C = -1.0;

        private static double x0 = 0;
        private static Vector<double> y0 = Vector<double>.Build.DenseOfArray(new double[] { 1.0, 1.0, A, 1.0 });

        private const double aCon = 0.0;
        private const double bCon = 5.0;


        private static Vector<double> Function(double x, Vector<double> y)
        {
            return Vector<double>.Build.DenseOfArray(new double[]
            {
                2 * x * Math.Pow(y[1], 1.0 / B) * y[3],
                2 * B * x * Math.Exp((B / C) * (y[2] - A)) * y[3],
                2 * C * x * y[3],
                -2 * x * Math.Log(y[0])
            });
        }

        public static Vector<double> RealFunction(double x)
        {
            return Vector<double>.Build.DenseOfArray(new double[]
            {
                Math.Exp(Math.Sin(x*x)),
                Math.Exp(B * Math.Sin(x*x)),
                C * Math.Sin(x*x) + A,
                Math.Cos(x * x)
            });
        }

        public static Vector<double> MyMethod(double h, double a = aCon, double b = bCon)
        {
            double c2 = 0.70;
            double a21 = c2;
            double b2 = 1.0 / (2.0 * c2);
            double b1 = 1.0 - 1.0 / (2.0 * c2);

            Vector<double> y = y0.Clone();
            int n = (int)Math.Ceiling((b - a) / h);

            double tmpX = x0;
            Vector<double> tmpY = y0.Clone();

            for (int i = 1; i <= n; ++i)
            {
                y = tmpY + h * (b1 * Function(tmpX, tmpY) + b2 * Function(tmpX + c2 * h, tmpY + h * a21 * Function(tmpX, tmpY)));
                tmpX = tmpX + h;
                tmpY = y;
            }
            return y;
        }

        public static Vector<double> HisMethod(double h, double a = aCon, double b = bCon)
        {
            Vector<double> y = Vector<double>.Build.Dense(y0.Count);
            int n = (int)Math.Ceiling((b - a) / h);

            double tmpX = x0;
            Vector<double> tmpY = y0.Clone();

            for (int i = 1; i <= n; ++i)
            {
                y = tmpY + h * Function(tmpX + h / 2.0, tmpY + (h / 2.0) * Function(tmpX, tmpY));
                tmpX = tmpX + h;
                tmpY = y;
            }
            return y;
        }

        public static double RungeRuleGlobal(double tol, double a = aCon, double b = bCon)
        {

            double h = Math.Pow(tol, 1.0 / 2.0);
            Vector<double> yn = MyMethod(h, a, b);
            Vector<double> y2n = MyMethod(h / 2.0, a, b);

            Vector<double> R2n = (y2n - yn) / (Math.Pow(2, 2) - 1.0);

            return h / 2.0 * Math.Pow(tol / R2n.L2Norm(), 1.0 / 2.0);
        }

        private static double InitialStep(double tol, double a = aCon, double b = bCon)
        {

            double h1 = Math.Pow(tol, 1.0 / 2.0);
            Vector<double> u = y0 + h1 * Function(x0, y0);
            double delta = Math.Pow(1.0 / Math.Max(a, b), 2.0 + 1.0) + Math.Pow(Function(x0 + h1, u).L2Norm(), 2.0 + 1.0);
            double h2 = Math.Pow(tol / delta, 1.0 / (2.0 + 1.0));
            return Math.Min(h1, h2);
        }
    }
}
