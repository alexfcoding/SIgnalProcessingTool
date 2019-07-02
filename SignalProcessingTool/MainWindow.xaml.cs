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
            double pi = Math.PI;
            double func = 0;
            
            double diameter = 0.038;
            double thickness = 0.003;
            double density = 7800;
            double speed = 0;
            double area = 0.3;
            double pipeLength = 10;
            double xCoordinate = pipeLength/1.6;
            
            double a1 = 10;
            double a2 = 10;
            double a3 = 5000000;
            double signalDuration = 1;
            double timeContact = 0.00004;
            double ti = 0.000022;
            double Fd = 44100;

            double modNumber = 500;

            double[] Cn = new double[500];
            double[] delta = new double[500];
            double[] modFrequency = new double[500];
            double[] Oi = new double[500];
            double[] sum = new double[88200];

            double E = 200 * Math.Pow(10, 9);
            double J = pi * Math.Pow(diameter, 3) * (double)thickness / 8;
            double Mass = density * pi * diameter * thickness;
            double a4 = E * J / Mass;
           
            int k = 0;

            for (int i = 1; i < modNumber; i = i + 2)
            {
                Cn[k] = pi / 2 * (2 * i + 1);
                delta[k] = 0.5 * (a1 * Math.Pow(i, 1) * Math.Pow(pi, 4) / Math.Pow(pipeLength, 4) + a2);
                modFrequency[k] = Math.Sqrt(a3 + a4 * Math.Pow(i, 4) * Math.Pow(pi, 4) / Math.Pow(pipeLength, 4));
                Oi[k] = Math.Sqrt(Math.Pow(modFrequency[k], 2) - Math.Pow(delta[k], 2));
                k++;
            }

            k = 0;

            int m = 0;

            for (int i = 1; i < signalDuration * Fd; i++)
            {
                for (int j = 1; j < modNumber; j = j + 2)
                {
                    sum[m] = sum[m] + 
                    Math.Pow(-1, (j - 1) / 2) * 
                    Math.Exp(-delta[k] * ti) * 
                    Math.Sin(j * pi * xCoordinate / pipeLength) /
                    (Math.Pow((Math.Pow(modFrequency[k], 2) - Math.Pow(pi, 2) / Math.Pow(timeContact, 2)), 2) + 
                    4 * Math.Pow(pi, 2) * Math.Pow(delta[k], 2) / Math.Pow(timeContact, 2)) *
                    ((2 * Math.Exp(delta[k] * timeContact) * 
                    delta[k] * pi / timeContact * Math.Cos(Oi[k] * timeContact) - 
                    pi * Math.Exp(delta[k] * timeContact) / timeContact / Oi[k] *
                    (2 * Math.Pow(delta[k], 2) - Math.Pow(modFrequency[k], 2) + 
                    Math.Pow(pi, 2) / timeContact / timeContact) * Math.Sin(Oi[k] * timeContact) +
                    2 * pi * delta[k] / timeContact) * Math.Cos(Oi[k] * ti) +
                    (pi * Math.Exp(delta[k] * timeContact) / timeContact / Oi[k] * (2 * Math.Pow(delta[k], 2) - 
                    Math.Pow(modFrequency[k], 2) + Math.Pow(pi, 2) / timeContact / timeContact) * Math.Cos(Oi[k] * timeContact) +
                    2 * Math.Exp(delta[k] * timeContact) * delta[k] * pi / timeContact * Math.Sin(Oi[k] * timeContact) + 
                    pi / timeContact / Oi[k] * (2 * Math.Pow(delta[k], 2) - Math.Pow(modFrequency[k], 2) +
                    Math.Pow(pi, 2) / timeContact / timeContact)) * Math.Sin(Oi[k] * ti));
                        
                    k++;                   
                }

                k = 0;
                m++;
                ti = ti + (double) (1 / Fd);
            }

            DrawImpulse(signalDuration * Fd - 1, sum);
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

    public class PointClass
    {
        public double xPoint { get; set; }
        public double yPoint { get; set; }
        public double yPointMax { get; set; }
    }
}
