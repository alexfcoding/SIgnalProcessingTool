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
using DSPLib;

namespace SignalProcessingTool
{
    using System.Collections.Generic;
    using OxyPlot;
    using System.Diagnostics;
    using OxyPlot.Series;
    using OxyPlot.Axes;
    using System.Collections.ObjectModel;
    using System.Numerics;

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
            Model1.Fd = 44100;
            Model1.ModNumber = 500;
            Model1.Cn = new double[500];
            Model1.Delta = new double[500];
            Model1.W = new double[500];
            Model1.Oi = new double[500];
            Model1.Sum = new double[44100];

            Model1.XCoordinate = Model1.PipeLength / 1.6;
            Model1.E = 200 * Math.Pow(10, 9);
            Model1.J = Model1.Pi * Math.Pow(Model1.Diameter, 3) * (double) Model1.Thickness / 8;
            Model1.Mass = Model1.Density * Model1.Pi * Model1.Diameter * Model1.Thickness;
            Model1.A4 = Model1.E * Model1.J / Model1.Mass;

            Model1.ComputeModel();

            Collection<PointClass> pointsImpulse = new Collection<PointClass>();
            DrawImpulse(pointsImpulse, Model1.Sum, false, Plot1);

            //double[] fftValues = FastFourier(Model1.Sum, false).Item1;
            //double[] fftFreq = FastFourier(Model1.Sum, false).Item2;
            double[] fftValues = FFTMathNumerics(Model1.Sum).Item1;
            double[] fftFreq = FFTMathNumerics(Model1.Sum).Item2;


            Collection<PointClass> pointsSpectrum = new Collection<PointClass>();
            // DrawImpulse(pointsSpectrum, fftValues, true, Plot2, fftFreq);
            DrawImpulse(pointsSpectrum, fftValues, false, Plot2);
        }

        void DrawImpulse(Collection<PointClass> pointCollection, double[] valuesArray, bool isSpectrum, OxyPlot.Wpf.Plot plotToDraw, double[]fftFreq = null)
        {
            plotToDraw.Series[0].ItemsSource = pointCollection;

            for (int i = 0; i < valuesArray.Length; i = i + 1)
            {
                if (isSpectrum == false)
                {
                    pointCollection.Add(
                    new PointClass
                    {
                        xPoint = i,
                        yPoint = valuesArray[i],
                    });
                }
                else
                {
                    pointCollection.Add(
                    new PointClass
                    {
                        xPoint = fftFreq[i],
                        yPoint = valuesArray[i],
                    });
                }
            }

            plotToDraw.InvalidatePlot(true);
        }

        Tuple<double[], double[]> FastFourier(double[] inputArray, bool logScale)
        {
            double[] tempArray = new double[2048];

            for (int i = 0; i < 2048; i++)
            {
                tempArray[i] = inputArray[i];
            }

            UInt32 length = 2048;
            double samplingRate = 44100;
            double[] spectrum = new double[4096];
            double[] wCoefs = DSP.Window.Coefficients(DSP.Window.Type.Hamming, length);
            double[] wInputData = DSP.Math.Multiply(tempArray, wCoefs);
            double wScaleFactor = DSP.Window.ScaleFactor.Signal(wCoefs);

            DSPLib.FFT fft = new DSPLib.FFT();

            fft.Initialize(length, length * 3);

            Complex[] cSpectrum = fft.Execute(wInputData);
            spectrum = DSP.ConvertComplex.ToMagnitude(cSpectrum);

            if (logScale == true)
            {
                spectrum = DSP.ConvertMagnitude.ToMagnitudeDBV(spectrum);

                for (int i = 0; i < spectrum.Length; i++)
                {
                    spectrum[i] -= 51;
                }
            }

            spectrum = DSP.Math.Multiply(spectrum, wScaleFactor);

            double[] freqSpan = fft.FrequencySpan(samplingRate);

            var tuple = new Tuple<double[], double[]>(spectrum, freqSpan);
            
            return tuple;
        }

        Tuple<double[], double[]> FFTMathNumerics (double[] inputArray)
        {
            double[] tempArray = new double[32768];

            for (int i = 0; i < 32768; i++)
            {
                tempArray[i] = inputArray[i];
            }

            Complex[] complexInput = new Complex[tempArray.Length];
            for (int i = 0; i < complexInput.Length; i++)
            {
                Complex tmp = new Complex(tempArray[i], 0);
                complexInput[i] = tmp;
            }

            MathNet.Numerics.IntegralTransforms.Fourier.Forward(complexInput);

            //do some stuff

            MathNet.Numerics.IntegralTransforms.Fourier.Inverse(complexInput);

            double[] outSamples = new double[complexInput.Length];

            for (int i = 0; i < outSamples.Length; i++)
                outSamples[i] = (double)complexInput[i].Real;

            double[] freqSpan = MathNet.Numerics.IntegralTransforms.Fourier.FrequencyScale(32768, 44100);

            var tuple = new Tuple<double[], double[]>(outSamples, freqSpan);

            return tuple;
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
            Fd = fd1;
            ModNumber = modNumber;
            Cn = cn1;
            Delta = delta;
            W = w1;
            Oi = oi1;
            Sum = sum;
            E = e1;
            J = j1;
            Mass = mass1;
            A4 = a4;
            Pi = pi;
        }

        public void ComputeModel()
        {
            int k = 0;

            for (int i = 1; i < ModNumber; i = i + 2)
            {
                Cn[k] = Pi / 2 * (2 * i + 1);
                Delta[k] = 0.5 * (A1 * Math.Pow(i, 1) * Math.Pow(Pi, 4) / Math.Pow(PipeLength, 4) + A2);
                W[k] = Math.Sqrt(A3 + A4 * Math.Pow(i, 4) * Math.Pow(Pi, 4) / Math.Pow(PipeLength, 4));
                Oi[k] = Math.Sqrt(Math.Pow(W[k], 2) - Math.Pow(Delta[k], 2));
                k++;
            }

            k = 0;

            double m1, m2, m3, rotator;

            for (int i = 1; i < SignalDuration * Fd; i++)
            {
                for (int j = 1; j < ModNumber; j = j + 2)
                {
                    rotator = Math.Pow(-1, (j - 1) / 2);

                    m1 = 2 * Math.Exp(Delta[k] * Tc) * Delta[k] * Pi / Tc;

                    m2 = Pi * Math.Exp(Delta[k] * Tc) / Tc / Oi[k] * (2 * Math.Pow(Delta[k], 2) - Math.Pow(W[k], 2) + Math.Pow(Pi, 2) / Math.Pow(Tc, 2));

                    m3 = Oi[k] * Tc;

                    Sum[i] = Sum[i] + rotator * Math.Exp(-Delta[k] * Ti) * Math.Sin(j * Pi * XCoordinate / PipeLength) /
                    (Math.Pow((Math.Pow(W[k], 2) - Math.Pow(Pi, 2) / Math.Pow(Tc, 2)), 2) + 4 * Math.Pow(Pi, 2) * Math.Pow(Delta[k], 2) / Math.Pow(Tc, 2)) *

                    ((m1 * Math.Cos(m3) - m2 * Math.Sin(m3) + 2 * Pi * Delta[k] / Tc) * Math.Cos(Oi[k] * Ti) + (m2 * Math.Cos(m3) + m1 * Math.Sin(m3) +
                    Pi / Tc / Oi[k] * (2 * Math.Pow(Delta[k], 2) - Math.Pow(W[k], 2) + Math.Pow(Pi, 2) / Math.Pow(Tc, 2))) * Math.Sin(Oi[k] * Ti));

                    k++;
                }

                k = 0;
                
                Ti = Ti + (double)(1 / Fd);
            }
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
        public double Fd { get; set; }
        public double ModNumber { get; set; }
        public double[] Cn { get; set; }
        public double[] Delta { get; set; }
        public double[] W { get; set; }
        public double[] Oi { get; set; }
        public double[] Sum { get; set; }
        public double E { get; set; }
        public double J { get; set; }
        public double Mass { get; set; }
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
