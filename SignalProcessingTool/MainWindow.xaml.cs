using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Forms.DataVisualization.Charting;
using System.Drawing;
using System.Diagnostics;


namespace SignalProcessingTool
{
    using System.Collections.Generic;
    using OxyPlot;
    using System.Diagnostics;
    using OxyPlot.Series;
    using OxyPlot.Axes;
    using System.Collections.ObjectModel;

    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void BtnCalcStart_Click(object sender, RoutedEventArgs e)
        {
            CalculateImpulse();
        }

        public void CalculateImpulse()
        {
            AcousticModel Model1 = new AcousticModel();

            Model1.Pi = Math.PI;
            Model1.Diameter = 0.038;
            Model1.Thickness = 0.003;
            Model1.Density = 7800;
            Model1.PipeLength = 10;
            
            Model1.A1 = 10;
            Model1.A2 = 10;
            Model1.A3 = 5000000;
            Model1.SignalDuration = 1;
            Model1.Tc = 0.00004;
            Model1.Ti = 0.000022;
            Model1.Fd1 = 44100;
            Model1.ModNumber = 500;
            Model1.Cn1 = new double[500];
            Model1.Delta = new double[500];
            Model1.W1 = new double[500];
            Model1.Oi1 = new double[500];
            Model1.Sum = new double[44100];

            Model1.XCoordinate = Model1.PipeLength / 1.6;
            Model1.E1 = 200 * Math.Pow(10, 9);
            Model1.J1 = Model1.Pi * Math.Pow(Model1.Diameter, 3) * (double) Model1.Thickness / 8;
            Model1.Mass1 = Model1.Density * Model1.Pi * Model1.Diameter * Model1.Thickness;
            Model1.A4 = Model1.E1 * Model1.J1 / Model1.Mass1;

            Model1.ComputeModel();

            DrawImpulse(Model1.SignalDuration * Model1.Fd1 - 1, Model1.Sum);
            //Plot2.Series[0].ItemsSource = Model1.Points;

        }

        void DrawImpulse(double xPointsCount, double[] valuesArray)
        {
            Collection<PointClass> Points = new Collection<PointClass>();

            Plot1.Series[0].ItemsSource = Points;

            for (int i = 0; i < xPointsCount; i = i + 1)
            {
                Points.Add(
                new PointClass
                {
                    xPoint = i,
                    yPoint = valuesArray[i],
                });
            }

            Plot1.InvalidatePlot(true);
        }


    }

    public class AcousticModel
    {
        public AcousticModel() { }

        public AcousticModel(double diameter, double thickness, double density, double pipeLength, double xCoordinate, double a1, double a2, double a3, double signalDuration, double tc, double ti, double fd1, double modNumber, double[] cn1, double[] delta, double[] w1, double[] oi1, double[] sum, double e1, double j1, double mass1, double a4, double pi)
        {
            Diameter = diameter;
            Thickness = thickness;
            Density = density;
            PipeLength = pipeLength;
            XCoordinate = xCoordinate;
            A1 = a1;
            A2 = a2;
            A3 = a3;
            SignalDuration = signalDuration;
            Tc = tc;
            Ti = ti;
            Fd1 = fd1;
            ModNumber = modNumber;
            Cn1 = cn1;
            Delta = delta;
            W1 = w1;
            Oi1 = oi1;
            Sum = sum;
            E1 = e1;
            J1 = j1;
            Mass1 = mass1;
            A4 = a4;
            Pi = pi;
        }

        public void ComputeModel()
        {
            int k = 0;

            for (int i = 1; i < ModNumber; i = i + 2)
            {
                Cn1[k] = Pi / 2 * (2 * i + 1);
                Delta[k] = 0.5 * (A1 * Math.Pow(i, 1) * Math.Pow(Pi, 4) / Math.Pow(PipeLength, 4) + A2);
                W1[k] = Math.Sqrt(A3 + A4 * Math.Pow(i, 4) * Math.Pow(Pi, 4) / Math.Pow(PipeLength, 4));
                Oi1[k] = Math.Sqrt(Math.Pow(W1[k], 2) - Math.Pow(Delta[k], 2));
                k++;
            }

            k = 0;

            double m1, m2, m3, rotator;

            for (int i = 1; i < SignalDuration * Fd1; i++)
            {
                for (int j = 1; j < ModNumber; j = j + 2)
                {
                    rotator = Math.Pow(-1, (j - 1) / 2);

                    m1 = 2 * Math.Exp(Delta[k] * Tc) * Delta[k] * Pi / Tc;

                    m2 = Pi * Math.Exp(Delta[k] * Tc) / Tc / Oi1[k] * (2 * Math.Pow(Delta[k], 2) - Math.Pow(W1[k], 2) + Math.Pow(Pi, 2) / Math.Pow(Tc, 2));

                    m3 = Oi1[k] * Tc;

                    Sum[i] = Sum[i] + rotator * Math.Exp(-Delta[k] * Ti) * Math.Sin(j * Pi * XCoordinate / PipeLength) /
                    (Math.Pow((Math.Pow(W1[k], 2) - Math.Pow(Pi, 2) / Math.Pow(Tc, 2)), 2) + 4 * Math.Pow(Pi, 2) * Math.Pow(Delta[k], 2) / Math.Pow(Tc, 2)) *

                    ((m1 * Math.Cos(m3) - m2 * Math.Sin(m3) + 2 * Pi * Delta[k] / Tc) * Math.Cos(Oi1[k] * Ti) + (m2 * Math.Cos(m3) + m1 * Math.Sin(m3) +
                    Pi / Tc / Oi1[k] * (2 * Math.Pow(Delta[k], 2) - Math.Pow(W1[k], 2) + Math.Pow(Pi, 2) / Math.Pow(Tc, 2))) * Math.Sin(Oi1[k] * Ti));

                    k++;
                }

                k = 0;
                
                Ti = Ti + (double)(1 / Fd1);
            }

            //DrawImpulse(signalDuration * Fd - 1, sum);
        }

        public double Diameter { get; set; }
        public double Thickness { get; set; }
        public double Density { get; set; }
        public double PipeLength { get; set; }
        public double XCoordinate { get; set; }
        public double A1 { get; set; }
        public double A2 { get; set; }
        public double A3 { get; set; }
        public double SignalDuration { get; set; }
        public double Tc { get; set; }
        public double Ti { get; set; }
        public double Fd1 { get; set; }
        public double ModNumber { get; set; }
        public double[] Cn1 { get; set; }
        public double[] Delta { get; set; }
        public double[] W1 { get; set; }
        public double[] Oi1 { get; set; }
        public double[] Sum { get; set; }
        public double E1 { get; set; }
        public double J1 { get; set; }
        public double Mass1 { get; set; }
        public double A4 { get; set; }
        public double Pi { get; set; }
    }

    public class PointClass
    {
        public double xPoint { get; set; }
        public double yPoint { get; set; }
        public double yPointMax { get; set; }
    }
}
