package NeuralNetwork
  annotation (uses(Modelica(version="2.2")), Documentation(info="<html>

This is the main package of the NeuralNetwork library.

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"));

  package BaseClasses
    "The basic element of the NeuralNetwork library is modeled"
    block NeuralNetworkLayer
      "This is the basic model for a neural network layer"

      parameter Integer numNeurons=1
        "It specifies the number of neurons which compose the laye"
           annotation(Dialog(group="Layer Data Definition"));
      parameter Integer numInputs=1
        "It specifies the number of inputs of the layer "
           annotation(Dialog(group="Layer Data Definition"));
      parameter Real weightTable[:,:] = [0,  0;  0,  0]
        "It is the weight table of the layer "                                             annotation(Dialog(group="Layer Data Definition"));
      parameter Real biasTable[:,:] = [0,  0]
        "It is the bias table of the layer"                                    annotation(Dialog(group="Layer Data Definition"));
      parameter NeuralNetwork.Types.ActivationFunction.Temp
        NeuronActivationFunction =  NeuralNetwork.Types.ActivationFunction.TanSig
        "It is the activation function of the layer"
      annotation(Dialog(group="Layer Data Definition"));

    protected
      extends Modelica.Blocks.Interfaces.MIMO(final nin=numInputs,final nout=numNeurons);

    equation
      if (NeuronActivationFunction == NeuralNetwork.Types.ActivationFunction.PureLin) then
        y = weightTable * u + biasTable[:,1];
      elseif (NeuronActivationFunction == NeuralNetwork.Types.ActivationFunction.TanSig) then
        y = Modelica.Math.tanh(weightTable * u + biasTable[:,1]);
      elseif (NeuronActivationFunction == NeuralNetwork.Types.ActivationFunction.LogSig) then
        y= NeuralNetwork.Utilities.LogSig(weightTable * u + biasTable[:,1]);
      elseif (NeuronActivationFunction == NeuralNetwork.Types.ActivationFunction.RadBas) then
        y= vector(NeuralNetwork.Utilities.RadBas(matrix(NeuralNetwork.Utilities.ElementWiseProduct(matrix(NeuralNetwork.Utilities.Dist(weightTable,matrix(u))),matrix(biasTable)))));
      end if;

      annotation (Icon(Bitmap(extent=[-100,94; 100,-98], name="Icons/NN-Symbol.png")),
                                Diagram(Text(
            extent=[-102,118; 102,62],
            style(color=3, rgbcolor={0,0,255}),
            string="NeuralNetworkLayer"),
                       Bitmap(extent=[-112,78; 114,-96], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[-100,100; 100,-100], style(
              color=3,
              rgbcolor={0,0,255},
              thickness=4))),
        DymolaStoredErrors,
        Documentation(info="<html>

<p>
This block models a Neural Network Layer.
</p>

<p>
A Neural Network Layer is specified by the following parameters
<UL>
<LI> <i>numNeurons</i>: it specifies the number of neurons which compose the layer (it is also equal to the rows numer of the weight and bias matrix and to the number of outputs of the layer;
<LI> <i>numInputs</i>: it specifies the number of inputs of the layer (it is also equal to the columns numer of the weight matrix;
<LI> <i>weightTable</i>: it is the weight table of the layer ([Number of Neurons x Number of Inputs]);
<LI> <i>biasTable</i>: it is the bias table of the layer ([Number of Neurons x 1]);
<LI> <i>NeuronActivationFunction</i>: it is the activation function of the layer; it can be equal to:
<UL>
<LI> NeuralNetwork.Types.ActivationFunction.TanSig = Hyperbolic tangent sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.LogSig = Logarithmic sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.RadBas = Radial basis activation function
<LI> NeuralNetwork.Types.ActivationFunction.PureLin = Linear activation function
</p>


<p>
To get the weight and bias table as modelica wants two different ways can be used:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder;
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>


</html>"));

    end NeuralNetworkLayer;

    annotation (Documentation(info="<html>

In this package the basic element of the NeuralNetwork library is modeled.

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"));
  end BaseClasses;

  package Networks
    "The complex elements of the NeuralNetwork library are modeled"
    block NeuralNetwork_TwoLayer "This block models a two layer Neural Network"

      BaseClasses.NeuralNetworkLayer neuralNetwork_HiddenLayer(
        numNeurons=HiddenLayer_numNeurons,
        numInputs=HiddenLayer_numInputs,
        weightTable=HiddenLayer_weightTable,
        biasTable=HiddenLayer_biasTable,
        NeuronActivationFunction=HiddenLayer_NeuronActivationFunction)
        annotation (extent=[-80,-24; -26,24]);
      BaseClasses.NeuralNetworkLayer neuralNetwork_OutputLayer(
        numNeurons=OutputLayer_numNeurons,
        numInputs=OutputLayer_numInputs,
        weightTable=OutputLayer_weightTable,
        biasTable=OutputLayer_biasTable,
        NeuronActivationFunction=OutputLayer_NeuronActivationFunction)
        annotation (extent=[20,-24; 72,24]);

      annotation (Diagram,
        Icon(
          Bitmap(extent=[-98,42; -4,-44], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[-98,40; -6,-38], style(color=3, rgbcolor={0,0,255})),
          Bitmap(extent=[6,42; 100,-44], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[6,40; 98,-38], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[-6,6; 6,0; -6,-6; -6,6], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0}))),
        Documentation(info="<html>

<p>
This block models a two layer Neural Network.
</p>

<p>
A two layer Neural Network is composed by two NeuralNetworkLayer (HiddenLayer_... and OutputLayer_...). Everyone is specified by the following parameters:
<UL>
<LI> <i>numNeurons</i>: it specifies the number of neurons which compose the layer (it is also equal to the rows numer of the weight and bias matrix and to the number of outputs of the layer;
<LI> <i>numInputs</i>: it specifies the number of inputs of the layer (it is also equal to the columns numer of the weight matrix;
<LI> <i>weightTable</i>: it is the weight table of the layer ([Number of Neurons x Number of Inputs]);
<LI> <i>biasTable</i>: it is the bias table of the layer ([Number of Neurons x 1]);
<LI> <i>NeuronActivationFunction</i>: it is the activation function of the layer; it can be equal to:
<UL>
<LI> NeuralNetwork.Types.ActivationFunction.TanSig = Hyperbolic tangent sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.LogSig = Logarithmic sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.RadBas = Radial basis activation function
<LI> NeuralNetwork.Types.ActivationFunction.PureLin = Linear activation function
</p>


<p>
To get the weight and bias table as modelica wants two different ways can be used:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder;
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>


</html>"));

      parameter Integer HiddenLayer_numNeurons=1
        "It specifies the number of neurons which compose the hidden layer";
      parameter Integer HiddenLayer_numInputs=1
        "It specifies the number of inputs of the hidden layer ";
      parameter Real HiddenLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the hidden layer ";
      parameter Real HiddenLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the hidden layer";
      parameter Types.ActivationFunction.Temp
        HiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the hidden layer";
      parameter Integer OutputLayer_numNeurons=1
        "It specifies the number of neurons which compose the output layer";
      parameter Integer OutputLayer_numInputs=1
        "It specifies the number of inputs of the output layer ";
      parameter Real OutputLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the output layer ";
      parameter Real OutputLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the output layer";
      parameter Types.ActivationFunction.Temp
        OutputLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the output layer";

    protected
        extends Modelica.Blocks.Interfaces.MIMO(final nin=neuralNetwork_HiddenLayer.numInputs,final nout
          =                                                                                              neuralNetwork_OutputLayer.numNeurons);

    equation
      connect(neuralNetwork_HiddenLayer.y, neuralNetwork_OutputLayer.u)
        annotation (points=[-23.3,0; 14.8,0],style(color=74, rgbcolor={0,0,127}));
      connect(u, neuralNetwork_HiddenLayer.u)
                                             annotation (points=[-120,0; -85.4,
            0], style(color=74, rgbcolor={0,0,127}));
      connect(neuralNetwork_OutputLayer.y, y) annotation (points=[74.6,0; 110,0],
          style(color=74, rgbcolor={0,0,127}));
    end NeuralNetwork_TwoLayer;

    block NeuralNetwork_ThreeLayer
      "This block models a three layer Neural Network"

      BaseClasses.NeuralNetworkLayer neuralNetworkLayer_FirstHiddenLayer(
        numNeurons=FirstHiddenLayer_numNeurons,
        numInputs=FirstHiddenLayer_numInputs,
        weightTable=FirstHiddenLayer_weightTable,
        biasTable=FirstHiddenLayer_biasTable,
        NeuronActivationFunction=FirstHiddenLayer_NeuronActivationFunction)
        annotation (extent=[-84,-14; -52,14]);
      BaseClasses.NeuralNetworkLayer neuralNetworkLayer_SecondHiddenLayer(
        numNeurons=SecondHiddenLayer_numNeurons,
        numInputs=SecondHiddenLayer_numInputs,
        weightTable=SecondHiddenLayer_weightTable,
        biasTable=SecondHiddenLayer_biasTable,
        NeuronActivationFunction=SecondHiddenLayer_NeuronActivationFunction)
        annotation (extent=[-22,-14; 12,14]);
      BaseClasses.NeuralNetworkLayer neuralNetworkLayer_OutputLayer(
        numNeurons=OutputLayer_numNeurons,
        numInputs=OutputLayer_numInputs,
        weightTable=OutputLayer_weightTable,
        biasTable=OutputLayer_biasTable,
        NeuronActivationFunction=OutputLayer_NeuronActivationFunction)
        annotation (extent=[50,-14; 80,14]);

      annotation (Diagram, Icon(
          Bitmap(extent=[-96,26; -32,-26], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[-92,24; -36,-24], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[-36,2; -28,0; -36,-2; -36,2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Bitmap(extent=[-32,24; 32,-28], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[-28,24; 28,-24], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[28,2; 36,0; 28,-2; 28,2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Bitmap(extent=[32,24; 96,-28], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[36,24; 92,-24], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[92,2; 100,0; 92,-2; 92,2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Polygon(points=[-100,2; -92,0; -100,-2; -100,2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1))),
        Documentation(info="<html>
<p>
This block models a three layer Neural Network.
</p>

<p>
A three layer Neural Network is composed by three NeuralNetworkLayer (FirstHiddenLayer_..., SecondHiddenLayer_... and OutputLayer_...). Everyone is specified by the following parameters:
<UL>
<LI> <i>numNeurons</i>: it specifies the number of neurons which compose the layer (it is also equal to the rows numer of the weight and bias matrix and to the number of outputs of the layer;
<LI> <i>numInputs</i>: it specifies the number of inputs of the layer (it is also equal to the columns numer of the weight matrix;
<LI> <i>weightTable</i>: it is the weight table of the layer ([Number of Neurons x Number of Inputs]);
<LI> <i>biasTable</i>: it is the bias table of the layer ([Number of Neurons x 1]);
<LI> <i>NeuronActivationFunction</i>: it is the activation function of the layer; it can be equal to:
<UL>
<LI> NeuralNetwork.Types.ActivationFunction.TanSig = Hyperbolic tangent sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.LogSig = Logarithmic sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.RadBas = Radial basis activation function
<LI> NeuralNetwork.Types.ActivationFunction.PureLin = Linear activation function
</p>


<p>
To get the weight and bias table as modelica wants two different ways can be used:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder;
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>


</html>"));

      parameter Integer FirstHiddenLayer_numNeurons=1
        "It specifies the number of neurons which compose the first hidden layer";
      parameter Integer FirstHiddenLayer_numInputs=1
        "It specifies the number of inputs of the first hidden layer ";
      parameter Real FirstHiddenLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the first hidden layer ";
      parameter Real FirstHiddenLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the first hidden layer";
      parameter Types.ActivationFunction.Temp
        FirstHiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the first hidden layer";
      parameter Integer SecondHiddenLayer_numNeurons=1
        "It specifies the number of neurons which compose the second hidden layer";
      parameter Integer SecondHiddenLayer_numInputs=1
        "It specifies the number of inputs of the second hidden layer ";
      parameter Real SecondHiddenLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the second hidden layer ";
      parameter Real SecondHiddenLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the second hidden layer";
      parameter Types.ActivationFunction.Temp
        SecondHiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the second hidden layer";
      parameter Integer OutputLayer_numNeurons=1
        "It specifies the number of neurons which compose the output layer";
      parameter Integer OutputLayer_numInputs=1
        "It specifies the number of inputs of the output layer ";
      parameter Real OutputLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the output layer ";
      parameter Real OutputLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the output layer";
      parameter Types.ActivationFunction.Temp
        OutputLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the output layer";

    protected
      extends Modelica.Blocks.Interfaces.MIMO(final nin=neuralNetwork_FirstHiddenLayer.numInputs,final nout
          =                                                                                                 neuralNetwork_OutputLayer.numNeurons);

    equation
      connect(u, neuralNetworkLayer_FirstHiddenLayer.u)
                                                  annotation (points=[-120,0;
            -92,0; -92,1.77636e-015; -87.2,1.77636e-015],
                    style(color=74, rgbcolor={0,0,127}));
      connect(neuralNetworkLayer_FirstHiddenLayer.y,
        neuralNetworkLayer_SecondHiddenLayer.u)
        annotation (points=[-50.4,1.77636e-015; -38,-3.1606e-021; -40,0; -25.4,
            1.77636e-015],                 style(color=74, rgbcolor={0,0,127}));
      connect(neuralNetworkLayer_SecondHiddenLayer.y, neuralNetworkLayer_OutputLayer.u)
        annotation (points=[13.7,1.77636e-015; 20.2,-3.1606e-021; 20,0; 47,
            1.77636e-015],              style(color=74, rgbcolor={0,0,127}));
      connect(neuralNetworkLayer_OutputLayer.y, y)
        annotation (points=[81.5,1.77636e-015; 84,1.77636e-015; 84,0; 110,0],
                                          style(color=74, rgbcolor={0,0,127}));
    end NeuralNetwork_ThreeLayer;

    block NeuralNetwork_RecurrentOneLayer
      "This block models a two layer Neural Network, with one recurrent layer"
      BaseClasses.NeuralNetworkLayer neuralNetwork_HiddenLayer(
        numNeurons=HiddenLayer_numNeurons,
        numInputs=HiddenLayer_numInputs+HiddenLayer_numNeurons,
        weightTable=HiddenLayer_weightTable,
        biasTable=HiddenLayer_biasTable,
        NeuronActivationFunction=HiddenLayer_NeuronActivationFunction)
        annotation (extent=[-28,-20; 14,20]);
      BaseClasses.NeuralNetworkLayer neuralNetwork_OutputLayer(
        numNeurons=OutputLayer_numNeurons,
        numInputs=OutputLayer_numInputs,
        weightTable=OutputLayer_weightTable,
        biasTable=OutputLayer_biasTable,
        NeuronActivationFunction=OutputLayer_NeuronActivationFunction)
        annotation (extent=[46,-20; 88,20]);
    protected
    extends Modelica.Blocks.Interfaces.MIMO(final nin=neuralNetwork_HiddenLayer.numInputs - neuralNetwork_HiddenLayer.numNeurons,final nout
          =                                                                                                    neuralNetwork_OutputLayer.numNeurons);
      Utilities.UnitDelayMIMO unitDelayMIMO(samplePeriod=samplePeriod,final nin=neuralNetwork_HiddenLayer.numNeurons, final nout = neuralNetwork_HiddenLayer.numNeurons)
        annotation (extent=[-52,48; -32,68], rotation=0);

      Modelica.Blocks.Routing.Multiplex2 multiplex2_1(final n1=neuralNetwork_HiddenLayer.numNeurons, final n2=neuralNetwork_HiddenLayer.numInputs - neuralNetwork_HiddenLayer.numNeurons)
        annotation (extent=[-80,-4; -56,18]);
      annotation (Diagram, Icon(
          Bitmap(extent=[-80,32; -16,-20], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[-76,32; -20,-16], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[-20,8; -12,6; -20,4; -20,8],  style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Bitmap(extent=[18,26; 82,-26], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[22,24; 78,-24], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[78,2; 86,0; 78,-2; 78,2],  style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Polygon(points=[-100,2; -92,0; -100,-2; -100,2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Line(points=[-92,0; -76,0; -92,0; -94,0], style(color=3, rgbcolor={0,0,
                  255})),
          Rectangle(extent=[-60,60; -38,44], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85})),
          Text(
            extent=[-56,56; -42,48],
            style(color=3, rgbcolor={0,0,255}),
            string="z^-1"),
          Line(points=[-12,6; 2,6], style(color=3, rgbcolor={0,0,255})),
          Line(points=[-6,6; -6,52; -38,52], style(color=3, rgbcolor={0,0,255})),
          Line(points=[-60,52; -88,52; -88,16; -76,16], style(color=3, rgbcolor={0,
                  0,255})),
          Line(points=[2,6; 22,6], style(color=3, rgbcolor={0,0,255})),
          Line(points=[86,0; 100,0], style(color=3, rgbcolor={0,0,255}))),
        Documentation(info="<html>

<p>
This block models a two layer Neural Network, with one recurrent layer.
</p>

<p>
A recurrent Neural Network is composed at least two NeuralNetworkLayer (HiddenLayer_... and OutputLayer_...). Everyone is specified by the following parameters:
<UL>
<LI> <i>numNeurons</i>: it specifies the number of neurons which compose the layer (it is also equal to the rows numer of the weight and bias matrix and to the number of outputs of the layer;
<LI> <i>numInputs</i>: it specifies the number of inputs of the layer (it is also equal to the columns numer of the weight matrix;
<LI> <i>weightTable</i>: it is the weight table of the layer ([Number of Neurons x Number of Inputs]);
<LI> <i>biasTable</i>: it is the bias table of the layer ([Number of Neurons x 1]);
<LI> <i>NeuronActivationFunction</i>: it is the activation function of the layer; it can be equal to:
<UL>
<LI> NeuralNetwork.Types.ActivationFunction.TanSig = Hyperbolic tangent sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.LogSig = Logarithmic sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.RadBas = Radial basis activation function
<LI> NeuralNetwork.Types.ActivationFunction.PureLin = Linear activation function
</UL>
</UL>
</p>

<p>
The network is called \"recurrent\" because usually the output of the hidden layers is used as an input of the same layer: obviously the output has to delayed. In this case, this is done using the block NeuralNetwork.Utilities.UnitDelayMIMO. The samplePeriod of this block is a parameter of the network.
</p>

<p>
The model is made so that the number of inputs to the recurrent layer has not to considered the recurrent inputs. For example, if the layer has 1 non-recurrent input and 5 neurons then the number of all inputs to the layer will be 6, but number of inputs which has to be inserted has to be equal to 1.
</p>

<p>
To get the weight and bias table as modelica wants two different ways can be used:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder;
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>


</html>"));

    public
      parameter Modelica.SIunits.Time samplePeriod=0.1
        "Sample period of the Recurrent Layer";
      parameter Integer HiddenLayer_numNeurons=1
        "It specifies the number of neurons which compose the hidden layer";
      parameter Integer HiddenLayer_numInputs=1
        "It specifies the number of inputs of the hidden layer";
      parameter Real HiddenLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the hidden layer ";
      parameter Real HiddenLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the hidden layer";
      parameter Types.ActivationFunction.Temp
        HiddenLayer_NeuronActivationFunction=
          NeuralNetwork.Types.ActivationFunction.TanSig
        "It is the activation function of the hidden layer";
      parameter Integer OutputLayer_numNeurons=1
        "It specifies the number of neurons which compose the output layer";
      parameter Integer OutputLayer_numInputs=1
        "It specifies the number of inputs of the output layer ";
      parameter Real OutputLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the output layer ";
      parameter Real OutputLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the output layer";
      parameter Types.ActivationFunction.Temp
        OutputLayer_NeuronActivationFunction=
          NeuralNetwork.Types.ActivationFunction.TanSig
        "It is the activation function of the output layer";

    equation
      connect(multiplex2_1.y, neuralNetwork_HiddenLayer.u) annotation (points=[-54.8,7;
            -54,6; -44,6; -44,0; -32.2,0],                                                                        style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(neuralNetwork_HiddenLayer.y, neuralNetwork_OutputLayer.u) annotation (points=[16.1,0;
            41.8,0],                                                                                           style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(u, multiplex2_1.u2) annotation (points=[-120,0; -104,0; -104,0.4;
            -82.4,0.4],                                               style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(neuralNetwork_OutputLayer.y, y) annotation (points=[90.1,0; 110,0],               style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));

      connect(neuralNetwork_HiddenLayer.y, unitDelayMIMO.u) annotation (points=[16.1,0;
            26,0; 26,40; -68,40; -68,58; -54,58],
                                             style(color=74, rgbcolor={0,0,127}));
      connect(unitDelayMIMO.y, multiplex2_1.u1) annotation (points=[-31,58; -18,
            58; -18,80; -96,80; -96,13.6; -82.4,13.6],
                                 style(color=74, rgbcolor={0,0,127}));
    end NeuralNetwork_RecurrentOneLayer;

    block NeuralNetwork_RecurrentTwoLayer
      "This block models a three layer Neural Network, with two recurrent layers"
      BaseClasses.NeuralNetworkLayer neuralNetwork_FirstHiddenLayer(
        numNeurons=FirstHiddenLayer_numNeurons,
        numInputs=FirstHiddenLayer_numInputs+FirstHiddenLayer_numNeurons,
        weightTable=FirstHiddenLayer_weightTable,
        biasTable=FirstHiddenLayer_biasTable,
        NeuronActivationFunction=FirstHiddenLayer_NeuronActivationFunction)
        annotation (extent=[-58,-4; -36,16]);

      annotation (Diagram, Icon(
          Bitmap(extent=[-86,28; -38,-14], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[-86,24; -40,-8],  style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[-40,4; -32,2; -40,0; -40,4],  style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Polygon(points=[-100,2; -92,0; -100,-2; -100,2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Rectangle(extent=[-74,54; -58,44], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85})),
          Text(
            extent=[-72,52; -60,46],
            style(color=3, rgbcolor={0,0,255}),
            string="z^-1"),
          Line(points=[-94,0; -86,0], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1)),
          Line(points=[-74,50; -92,50; -92,16; -86,16], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1)),
          Rectangle(extent=[-4,-36; 12,-46], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85})),
          Text(
            extent=[-2,-38; 10,-44],
            style(color=3, rgbcolor={0,0,255}),
            string="z^-1"),
          Bitmap(extent=[-22,18; 26,-24], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[-22,14; 24,-18],  style(color=3, rgbcolor={0,0,255})),
          Line(points=[-32,2; -22,2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1)),
          Line(points=[-28,2; -28,50; -58,50], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1)),
          Line(points=[-4,-40; -30,-40; -30,-8; -22,-8], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1)),
          Bitmap(extent=[42,20; 90,-22], name="Icons/NN-Symbol.png"),
          Rectangle(extent=[42,14; 88,-18],   style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[24,0; 32,-2; 24,-4; 24,0],       style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Polygon(points=[88,0; 96,-2; 88,-4; 88,0],       style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0},
              fillPattern=1)),
          Line(points=[32,-2; 42,-2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1)),
          Line(points=[96,-2; 100,-2], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1)),
          Line(points=[34,-2; 34,-40; 12,-40], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=47,
              rgbfillColor={255,170,85},
              fillPattern=1))),
        Documentation(info="<html>

<p>
This block models a three layer Neural Network, with two recurrent layers.
</p>

<p>
A recurrent Neural Network is composed at least two NeuralNetworkLayer (HiddenLayer_... and OutputLayer_...). Everyone is specified by the following parameters:
<UL>
<LI> <i>numNeurons</i>: it specifies the number of neurons which compose the layer (it is also equal to the rows numer of the weight and bias matrix and to the number of outputs of the layer;
<LI> <i>numInputs</i>: it specifies the number of inputs of the layer (it is also equal to the columns numer of the weight matrix;
<LI> <i>weightTable</i>: it is the weight table of the layer ([Number of Neurons x Number of Inputs]);
<LI> <i>biasTable</i>: it is the bias table of the layer ([Number of Neurons x 1]);
<LI> <i>NeuronActivationFunction</i>: it is the activation function of the layer; it can be equal to:
<UL>
<LI> NeuralNetwork.Types.ActivationFunction.TanSig = Hyperbolic tangent sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.LogSig = Logarithmic sigmoid activation function
<LI> NeuralNetwork.Types.ActivationFunction.RadBas = Radial basis activation function
<LI> NeuralNetwork.Types.ActivationFunction.PureLin = Linear activation function
</UL>
</UL>
</p>

<p>
The network is called \"recurrent\" because usually the output of the hidden layers is used as an input of the same layer: obviously the output has to delayed. In this case, this is done using the block NeuralNetwork.Utilities.UnitDelayMIMO. The samplePeriod of this block is a parameter of the network.
</p>

<p>
The model is made so that the number of inputs to the recurrent layer has not to considered the recurrent inputs. For example, if the layer has 1 non-recurrent input and 5 neurons then the number of all inputs to the layer will be 6, but number of inputs which has to be inserted has to be equal to 1.
</p>

<p>
To get the weight and bias table as modelica wants two different ways can be used:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder;
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>


</html>"));
      BaseClasses.NeuralNetworkLayer neuralNetwork_OutputLayer(
        numNeurons=OutputLayer_numNeurons,
        numInputs=OutputLayer_numInputs,
        weightTable=OutputLayer_weightTable,
        biasTable=OutputLayer_biasTable,
        NeuronActivationFunction=OutputLayer_NeuronActivationFunction)
        annotation (extent=[70,-10; 90,10]);
      BaseClasses.NeuralNetworkLayer neuralNetwork_SecondHiddenLayer(
        numNeurons=SecondHiddenLayer_numNeurons,
        numInputs=SecondHiddenLayer_numInputs+SecondHiddenLayer_numNeurons,
        weightTable=SecondHiddenLayer_weightTable,
        biasTable=SecondHiddenLayer_biasTable,
        NeuronActivationFunction=SecondHiddenLayer_NeuronActivationFunction)
        annotation (extent=[18,-10; 40,10]);

    protected
      extends Modelica.Blocks.Interfaces.MIMO(nin=neuralNetwork_FirstHiddenLayer.numInputs-neuralNetwork_FirstHiddenLayer.numNeurons,nout=neuralNetwork_OutputLayer.numNeurons);

      Modelica.Blocks.Routing.Multiplex2 multiplex2_1(final n1= neuralNetwork_FirstHiddenLayer.numNeurons, final n2= neuralNetwork_FirstHiddenLayer.numInputs-neuralNetwork_FirstHiddenLayer.numNeurons)
        annotation (extent=[-88,-4; -68,16]);

      Utilities.UnitDelayMIMO unitDelayMIMO(
        samplePeriod=samplePeriod,
        final nin=neuralNetwork_FirstHiddenLayer.numNeurons,
        final nout =                                                                                           neuralNetwork_FirstHiddenLayer.numNeurons)
        annotation (extent=[-76,44; -56,64], rotation=180);
      Utilities.UnitDelayMIMO unitDelayMIMO1(samplePeriod=samplePeriod, final nin=neuralNetwork_SecondHiddenLayer.numNeurons,final nout
          =                                                                                                    neuralNetwork_SecondHiddenLayer.numNeurons)
        annotation (extent=[4,-50; 24,-30], rotation=180);

      Modelica.Blocks.Routing.Multiplex2 multiplex2_2(final n1=
            neuralNetwork_FirstHiddenLayer.numInputs -
            neuralNetwork_FirstHiddenLayer.numNeurons, final n2=
            neuralNetwork_FirstHiddenLayer.numInputs -
            neuralNetwork_SecondHiddenLayer.numNeurons)
        annotation (extent=[-12,-10; 8,10]);

    public
      parameter Modelica.SIunits.Time samplePeriod=0.1
        "Sample period of the Recurrent Layer";
      parameter Integer FirstHiddenLayer_numNeurons=1
        "It specifies the number of neurons which compose the first hidden layer";
      parameter Integer FirstHiddenLayer_numInputs=1
        "It specifies the number of inputs of the first hidden layer ";
      parameter Real FirstHiddenLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the first hidden layer ";
      parameter Real FirstHiddenLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the first hidden layer";
      parameter Types.ActivationFunction.Temp
        FirstHiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the first hidden layer";
      parameter Integer SecondHiddenLayer_numNeurons=1
        "It specifies the number of neurons which compose the second hidden layer";
      parameter Integer SecondHiddenLayer_numInputs=1
        "It specifies the number of inputs of the second hidden layer ";
      parameter Real SecondHiddenLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the second hidden layer ";
      parameter Real SecondHiddenLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the second hidden layer";
      parameter Types.ActivationFunction.Temp
        SecondHiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the second hidden layer";
      parameter Integer OutputLayer_numNeurons=1
        "It specifies the number of neurons which compose the output layer";
      parameter Integer OutputLayer_numInputs=1
        "It specifies the number of inputs of the output layer ";
      parameter Real OutputLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the output layer ";
      parameter Real OutputLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the output layer";
      parameter Types.ActivationFunction.Temp
        OutputLayer_NeuronActivationFunction=NeuralNetwork.Types.
          ActivationFunction.TanSig
        "It is the activation function of the output layer";

    equation
      connect(multiplex2_1.y, neuralNetwork_FirstHiddenLayer.u) annotation (points=[-67,6;
            -60.2,6],                                                                                  style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(u, multiplex2_1.u2) annotation (points=[-120,0; -90,0],                 style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(neuralNetwork_OutputLayer.y, y) annotation (points=[91,0; 110,0], style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(multiplex2_2.y, neuralNetwork_SecondHiddenLayer.u) annotation (points=[9,0;
            15.8,0],                                                                                  style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(neuralNetwork_SecondHiddenLayer.y, neuralNetwork_OutputLayer.u)
        annotation (points=[41.1,0; 68,0],                 style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=0,
          rgbfillColor={0,0,0},
          fillPattern=1));
      connect(neuralNetwork_FirstHiddenLayer.y, multiplex2_2.u1) annotation (points=
           [-34.9,6; -14,6], style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=47,
          rgbfillColor={255,170,85},
          fillPattern=1));
      connect(neuralNetwork_FirstHiddenLayer.y, unitDelayMIMO.u) annotation (
          points=[-34.9,6; -24,6; -24,54; -54,54], style(color=74, rgbcolor={0,
              0,127}));
      connect(unitDelayMIMO.y, multiplex2_1.u1) annotation (points=[-77,54; -96,
            54; -96,12; -90,12], style(color=74, rgbcolor={0,0,127}));
      connect(neuralNetwork_SecondHiddenLayer.y, unitDelayMIMO1.u) annotation (
          points=[41.1,0; 54,0; 54,-40; 26,-40], style(color=74, rgbcolor={0,0,
              127}));
      connect(unitDelayMIMO1.y, multiplex2_2.u2) annotation (points=[3,-40; -22,
            -40; -22,-6; -14,-6], style(color=74, rgbcolor={0,0,127}));
    end NeuralNetwork_RecurrentTwoLayer;

    block NeuralNetwork_RadialBasis
      "This block models a Radial Basis Neural Network"

      BaseClasses.NeuralNetworkLayer RadialBasis_HiddenLayer(
        numNeurons=HiddenLayer_numNeurons,
        numInputs=HiddenLayer_numInputs,
        weightTable=HiddenLayer_weightTable,
        biasTable=HiddenLayer_biasTable,
        NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.RadBas)
        annotation (extent=[-80,-24; -26,24]);
      BaseClasses.NeuralNetworkLayer PureLin_OutputLayer(
        numNeurons=OutputLayer_numNeurons,
        numInputs=OutputLayer_numInputs,
        weightTable=OutputLayer_weightTable,
        biasTable=OutputLayer_biasTable,
        NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.PureLin)
        annotation (extent=[20,-24; 72,24]);

      annotation (Diagram(
          Bitmap(extent=[4,22; 86,-26], name="Icons/NNlinear-Symbol.png"),
          Rectangle(extent=[20,24; 72,-24], style(color=0, rgbcolor={0,0,0}))),
        Icon(
          Rectangle(extent=[-98,40; -6,-38], style(color=3, rgbcolor={0,0,255})),
          Rectangle(extent=[6,40; 98,-38], style(color=3, rgbcolor={0,0,255})),
          Polygon(points=[-6,6; 6,0; -6,-6; -6,6], style(
              color=3,
              rgbcolor={0,0,255},
              fillColor=0,
              rgbfillColor={0,0,0})),
          Bitmap(extent=[-98,38; -6,-36], name="Icons/NNradial-Symbol.png"),
          Bitmap(extent=[8,38; 100,-38], name="Icons/NNlinear-Symbol.png")),
        Documentation(info="<html>

<p>
This block models a Radial Basis Neural Network.
</p>

<p>
A Radial Basis Neural Network is composed by two NeuralNetworkLayer (HiddenLayer_... and OutputLayer_...). Everyone is specified by the following parameters:
<UL>
<LI> <i>numNeurons</i>: it specifies the number of neurons which compose the layer (it is also equal to the rows numer of the weight and bias matrix and to the number of outputs of the layer;
<LI> <i>numInputs</i>: it specifies the number of inputs of the layer (it is also equal to the columns numer of the weight matrix;
<LI> <i>weightTable</i>: it is the weight table of the layer ([Number of Neurons x Number of Inputs]);
<LI> <i>biasTable</i>: it is the bias table of the layer ([Number of Neurons x 1]);
</UL>
The activation function of the layers is fixed: the HiddenLayer uses a RadBas anctivation function and the OutputLayer uses a PureLin anctivation function.
</p>


<p>
To get the weight and bias table as modelica wants two different ways was used:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder;
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>


</html>"));

      parameter Integer HiddenLayer_numNeurons=1
        "It specifies the number of neurons which compose the hidden layer";
      parameter Integer HiddenLayer_numInputs=1
        "It specifies the number of inputs of the hidden layer ";
      parameter Real HiddenLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the hidden layer ";
      parameter Real HiddenLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the hidden layer";
      parameter Integer OutputLayer_numNeurons=1
        "It specifies the number of neurons which compose the output layer";
      parameter Integer OutputLayer_numInputs=1
        "It specifies the number of inputs of the output layer ";
      parameter Real OutputLayer_weightTable[:,:]=[0,0; 0,0]
        "It is the weight table of the output layer ";
      parameter Real OutputLayer_biasTable[:,:]=[0,0]
        "It is the bias table of the output layer";

    protected
        extends Modelica.Blocks.Interfaces.MIMO(final nin=RadialBasis_HiddenLayer.numInputs,final nout = PureLin_OutputLayer.numNeurons);

    equation
      connect(RadialBasis_HiddenLayer.y, PureLin_OutputLayer.u)
        annotation (points=[-23.3,0; 14.8,0],style(color=74, rgbcolor={0,0,127}));
      connect(u, RadialBasis_HiddenLayer.u)  annotation (points=[-120,0; -85.4,
            0], style(color=74, rgbcolor={0,0,127}));
      connect(PureLin_OutputLayer.y, y)       annotation (points=[74.6,0; 110,0],
          style(color=74, rgbcolor={0,0,127}));
    end NeuralNetwork_RadialBasis;
    annotation (Documentation(info="<html>

In this package the complex elements of the NeuralNetwork library, based on the NeuralNetwork.BaseClasses.NeuralNetworkLayer are modeled.

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"));
  end Networks;

  package Utilities
    "some functions and models which are indispensable to model and to use some elements of the NeuralNetwork library"
    function LogSig "Logarithmic sigmoid activation function"
      input Real u "Input of the function";
      output Real y "Output of the function";
    algorithm
      y :=1/(1 + Modelica.Math.exp(-u));
      annotation (Documentation(info="<html>
<p>
This function calculate the logaritmic sigmoid transfer function.
</p>
<p>
The algorithm used when the input is n is: logsig(n) = a = 1 / (1 + exp(-n)).
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</html>"));
    end LogSig;

    function Dist
      "Calculate the euclideal distance between two matrix of real numbers: WeightPoint(AxB) , InputPoint(BxC) -> Distance(AxC)"
      input Real WeightPoint[:,:] "It's the matrix of the weight points";
      input Real InputPoint[:,:] "It's the matrix of the input points";
      output Real Distance[size(WeightPoint,1),size(InputPoint,2)]
        "It's the matrix of the distances";

    protected
        Real temp;

    algorithm
      for i in 1:size(WeightPoint,1) loop
        for j in 1:size(InputPoint,2) loop
          temp:=0;
          for k in 1:size(WeightPoint,2) loop
            temp:=(WeightPoint[i,k]-InputPoint[k,j])^2 + temp;
          end for;
          Distance[i,j]:=sqrt(temp);
        end for;
      end for;

      annotation (Documentation(info="<html>
<p>
This function calculate the euclideal distance between two matrix of real numbers. This function is used to calculate the activation function of a radial basis neuron: the first input is the WeightPoint(AxB) and the second is the InputPoint(BxC); the output is the Distance(AxC). </p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</html>"));
    end Dist;

    function RadBas "Radial Basis activation function"
      input Real[:,:] n "Input of the function";
      output Real[size(n,1),size(n,2)] a "Output of the function";

    algorithm
      a := exp(-(NeuralNetwork.Utilities.ElementWiseProduct(n,n)));
      annotation (Documentation(info="<html>
<p>
This function calculate the radial basis transfer function.
</p>
<p>
The algorithm used when the input is n is: radbas(n) = a = exp(-n^2).
</p>
<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</html>"));
    end RadBas;

    function ElementWiseProduct "ElementWiseProduct of two matrix"
      input Real[:,:] u1 "First input matrix";
      input Real[:,:] u2 "Second input matrix";
      output Real[size(u1,1),size(u1,2)] y "Output matrix";

    algorithm
      for i in 1:size(u1,1) loop
        for j in 1:size(u1,2) loop
          y[i,j] := u1[i,j]*u2[i,j];
        end for;
      end for;

      annotation (Documentation(info="<html>
<p>
This function calculate the element wise product of two matrix.
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</html>"));
    end ElementWiseProduct;

    block SamplerMIMO "Ideal sampling of continuous MIMO signals"
      extends Modelica.Blocks.Interfaces.DiscreteMIMO;

      annotation (
        Icon(
          Ellipse(extent=[-25, -10; -45, 10], style(color=3, fillColor=7)),
          Ellipse(extent=[45, -10; 25, 10], style(color=73, fillColor=7)),
          Line(points=[-100, 0; -45, 0], style(color=3)),
          Line(points=[45, 0; 100, 0], style(color=73)),
          Line(points=[-35, 0; 30, 35], style(color=3))),
        Diagram(
          Ellipse(extent=[-25, -10; -45, 10], style(color=3, fillColor=7)),
          Ellipse(extent=[45, -10; 25, 10], style(color=73, fillColor=7)),
          Line(points=[-100, 0; -45, 0], style(color=3)),
          Line(points=[45, 0; 100, 0], style(color=73)),
          Line(points=[-35, 0; 30, 35], style(color=3))),
        Documentation(info="<HTML>
<p>
Samples the continues input signal with a sampling rate defined
via parameter <b>samplePeriod</b>. The input and output are vectors.
</p>
</HTML>
"));
    equation
      when {sampleTrigger, initial()} then
        y = u;
      end when;
    end SamplerMIMO;

    block UnitDelayMIMO "Unit Delay Block for MIMO signals"
      extends Modelica.Blocks.Interfaces.DiscreteMIMO;
      parameter Real y_start[nin]=fill(0.0,nin)
        "Initial value of output signal";

      annotation (
        Coordsys(
          extent=[-100, -100; 100, 100],
          grid=[2, 2],
          component=[20, 20]),
        Window(
          x=0.24,
          y=0.09,
          width=0.6,
          height=0.6),
        Documentation(info="<html>
<p>
This block describes a unit delay:
</p>
<pre>
          1
     y = --- * u
          z
</pre>
<p>
that is, the output signal y is the input signal u of the
previous sample instant. Before the second sample instant,
the output y is identical to parameter yStart. The variable y and u are vectors.
</p>
<p><b>Release Notes:</b></p>
<ul>
<li><i>June 13, 2000</i>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized, based on a model of Dieter Moormann.
</li>
</ul>
</HTML>
"),     Icon(
          Line(points=[-30, 0; 30, 0]),
          Text(extent=[-90, 10; 90, 90], string="1"),
          Text(extent=[-90, -10; 90, -90], string="z")),
        Diagram(
          Rectangle(extent=[-60, 60; 60, -60], style(color=73)),
          Text(extent=[-160, 10; -140, -10], string="u"),
          Text(extent=[115, 10; 135, -10], string="y"),
          Line(points=[-100, 0; -60, 0], style(color=73)),
          Line(points=[60, 0; 100, 0], style(color=73)),
          Line(points=[40, 0; -40, 0], style(color=0)),
          Text(
            extent=[-55, 55; 55, 5],
            string="1",
            style(color=0)),
          Text(
            extent=[-55, -5; 55, -55],
            string="z",
            style(color=0))),
        DymolaStoredErrors);
    equation
      when sampleTrigger then
        y = pre(u);
      end when;

    initial equation
        y = y_start;
    end UnitDelayMIMO;
    annotation (Documentation(info="<html>

In this package there are some functions and models which are indispensable to model and to use some elements of the NeuralNetwork library.

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"));
  end Utilities;

  package Types "Constants and types with choices, especially to build menus"
    extends Modelica.Icons.Library;

    package ActivationFunction
      "Type, constants and menu choices to define the activation function of a neural network layer"
      annotation (preferedView="text");
      extends Modelica.Icons.Library;

      constant Integer PureLin=0;
      constant Integer TanSig=1;
      constant Integer LogSig=2;
      constant Integer RadBas=3;

      type Temp
        "Temporary type of ActivationFunction with choices for menus (until enumerations are available)"

        extends Integer;
        annotation (choices(
          choice=NeuralNetwork.Types.ActivationFunction.PureLin
              "PureLin activation function",
          choice=NeuralNetwork.Types.ActivationFunction.TanSig
              "TanSig activation function",
          choice=NeuralNetwork.Types.ActivationFunction.LogSig
              "LogSig activation function",
          choice=NeuralNetwork.Types.ActivationFunction.RadBas
              "RadialBasis activation function"));
      end Temp;
    end ActivationFunction;

    annotation (Documentation(info="<HTML>
<p>
In this package <b>types</b> are defined that are used
in library NeuralNetwork. The types have additional annotation choices
definitions that define the menus to be built up in the graphical
user interface when the type is used as parameter in a declaration.
</p>
</HTML>"));
  end Types;

  package Examples "Some examples of using the NeuralNetwork library"



    model RecurrentNeuralNetwork "An example of a Recurrent NeuralNetwork"
      annotation (Diagram(
          Rectangle(extent=[-96,98; 84,54], style(color=3, rgbcolor={0,0,255})),
          Rectangle(extent=[-96,50; 84,-2], style(color=3, rgbcolor={0,0,255})),
          Text(
            extent=[50,66; 82,50],
            style(color=3, rgbcolor={0,0,255}),
            string="extractData.m"),
          Text(
            extent=[52,6; 78,0],
            style(color=3, rgbcolor={0,0,255}),
            string="DataFiles")),
          Documentation(info="<html>
<p>



This model show an example of a Recurrent NeuralNetwork created with the <i>Neural Network library</i>.
A Recurrent Neural Network is composed by a layer of neurons which adopts the TanSig activation function and a layer of newrons which perform a LogSig transformation.
</p>

<p>
This example is made using the following MatLab commands:
<UL>
<LI> P = round(rand(1,80));
<LI> T = [0 (P(1:end-1)+P(2:end) == 2)];
<LI> Pseq = con2seq(P);
<LI> Tseq = con2seq(T);
<LI> net = newelm([0 1],[20 1],{'tansig','logsig'});
<LI> net.trainParam.epochs = 500;
<LI> net = train(net,Pseq,Tseq);
<LI> Y = sim(net,Pseq);
<LI> z = seq2con(Y);
<LI> figure
<LI> plot(T)
<LI> hold
<LI> plot(z{1,1},'r')
<LI> time=(1:length(P))*0.01;
<LI> in=[time' P'];
<LI> out=[time' z{1,1}'];
<LI> tanL_weight=[net.LW{1} net.IW{1,1}];
<LI> tanL_bias=net.b{1};
<LI> logL_weight=net.LW{2,1};
<LI> logL_bias=net.b{2};
<LI> typing this command in matlab: save testData_RecurrentNN -V4 in out tanL_weight tanL_bias logL_weight logL_bias
</UL>

These commands create and train a new Recurrent Neural Network (net). The detail about the network created can be see exploring the network object created by matlab. Which is important for our example are these informations:
<UL>
<LI> The first layer of the network, the TanSig one, expect 1 input (size of net.IW{1});
<LI> The first layer of the network, the TanSig one, is made by 20 newrons (size of net.IW{1});
<LI> The weight table of the first layer is net.IW{1};
<LI> The bias table of the first layer is net.b{1};
<LI> The second layer of the network, the LogSig one, expect 20 inputs (size of net.LW{2});
<LI> The second layer of the network, the LogSig one, is made by 1 newron (size of net.LW{2});
<LI> The weight table of the second layer is net.LW{2};
<LI> The bias table of the second layer is net.b{2};</UL>
</p>

<p>
In this example the parameters of the network are specified using two different way:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder
<UL>
<LI> extractData('dataHidden.txt','HiddenLayer',[net.LW{1},net.IW{1,1}],net.b{1},'tan',1)
<LI> extractData('dataOutput.txt','OutputLayer',[net.LW{2,1}],net.b{2},'log')
</UL>
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p>
The model is simulated on the same data used in MatLab so it is possible to test if the model gives the expected output.
<p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"),
        DymolaStoredErrors);

      Modelica.Blocks.Sources.CombiTimeTable in1(
        tableOnFile=true,
        tableName="in",
        fileName="testData_RecurrentNN.mat")
                              annotation (extent=[-92,66; -72,86]);
      Networks.NeuralNetwork_RecurrentOneLayer neuralNetwork_RecurrentOneLayer(
    HiddenLayer_numNeurons= 20,
    HiddenLayer_numInputs= 1,
    HiddenLayer_weightTable = [0.26275001, -0.16191802, 0.33544231, -0.66194073, 0.21904685, -0.16924848, 0.15118808, 0.45364439, 0.34106662, -0.31665220, -0.33530923, -0.00830987, -0.49713462, -0.29487471, 0.45460882, 0.09279143, -0.41925201, 0.51585397, 0.26022429, 0.22105398, 0.07245009;
             0.47748724, 0.24471344, 0.04418866, -0.31418557, -0.43608327, -0.39776726, -0.29259645, -0.47107238, -0.46854408, 0.18863103, 0.09062813, -0.16764313, 0.13236413, 0.29624859, 0.06588751, -0.46729562, 0.56444916, -0.13898466, 0.22884767, 0.08732982, 1.57082665;
             -0.53946413, -0.24624713, -0.44311691, 0.29359242, 0.41486128, -0.42766706, -0.60723227, -0.07691763, -0.42155191, 0.19495585, 0.67885364, -0.32813799, -0.15680672, 0.10108443, 0.14432761, -0.37424897, -0.14723787, 0.21596783, 0.42272350, 0.45358671, 0.75448817;
             -0.43558988, 0.20801745, 0.47580266, 0.04886076, 0.03132738, 0.51826319, -0.08065257, 0.18225986, 0.42770191, 0.15098632, -0.17448438, -0.43156093, 0.45415346, -0.15528803, 0.38341289, 0.39940543, -0.51180571, -0.41594198, -0.18700632, -0.55705675, 0.60134343;
             0.51212369, -0.31555185, -0.41986506, -0.11083167, 0.51118392, -0.24150784, 0.19873041, 0.01947013, -0.44737550, 0.06134588, 0.44704796, 0.43971255, 0.32949160, 0.36189423, -0.51019270, -0.54672246, 0.16118876, 0.39814575, -0.00794786, 0.00643280, -0.55373787;
             0.26566717, 0.23048444, 0.55000321, -0.32476688, 0.26738670, -0.11340111, -0.66324960, -0.09368310, -0.73587549, 0.57716097, -0.13017588, 0.42409698, -0.79985995, 0.44380129, 0.50804785, 0.49361725, -0.22327209, 0.10077041, 0.90127056, -0.01044839, 0.43941547;
             -0.04541634, 0.29215476, -0.04054864, 0.18681246, 0.01061152, 0.44960962, -0.13599946, -0.09976167, 0.63956535, 0.35527306, 0.11204062, -0.49098218, 0.43469135, 0.28377814, -0.62456669, -0.60588339, 1.01100589, -0.37064655, -0.24470628, -0.04851647, 2.69864857;
             0.27367867, 0.28826794, 0.56599901, 0.28778309, 0.67067080, -0.00904424, 0.01121982, 0.49255823, 0.12831361, -0.19427086, -0.33949739, -0.58390696, -0.05498601, -0.12386989, -0.06521530, 0.00257663, 0.27749025, -0.23860212, -0.03878670, 0.53882530, -1.45076504;
             0.36989121, 0.16797370, -0.25615667, 0.49104078, -0.46136375, -0.50742623, -0.74435841, 0.67051434, -0.31622099, 0.37964829, -0.17534188, 0.04632838, 0.01211665, 0.37458385, 0.12905153, 0.33669597, -0.41358548, 0.36516050, 0.77577175, -0.21113382, 0.76711422;
             -0.07437908, 0.81211820, 0.29548513, 0.52706713, -0.71856517, 0.06780260, 0.30583533, -0.05900166, 0.27426431, 0.00789060, -0.10847559, 0.30830462, -0.00314673, -0.58629113, 0.29402894, 0.28767967, -0.12703552, -0.44551804, -0.37831514, -0.00970070, 0.10760177;
             0.21721939, -0.28026271, -0.35586056, 0.61168539, 0.17818255, 0.45868637, -0.50201209, -0.06075839, 0.45277062, 0.17384786, 0.13855737, -0.29755874, -0.39230605, 0.33367336, 0.50887936, 0.08526675, -0.10849225, -0.35230099, 0.07288201, 0.57777087, -1.02923206;
             0.18831270, -0.30223310, 0.16282453, 0.33876719, -0.13049841, 0.00403408, -0.54464292, 0.74512772, -0.40495147, -0.42598938, -0.22788284, -0.08144722, 0.02981773, 0.18243560, -0.27733563, -0.18002333, -0.23266159, 0.31085832, 0.23095722, -0.54703175, -1.25447005;
             0.51688792, 0.34488904, 0.09387091, 0.05615853, -0.14497822, 0.72810316, 0.43556794, -0.02499796, -0.24701200, -0.20505213, 0.48974231, -0.22499195, -0.31439278, 0.44166774, 0.28113559, 0.15622612, -0.38469002, -0.31565089, -0.00824059, -0.14469364, 1.19655615;
             -0.54495404, -0.55356573, -0.25137476, 0.47597406, -0.19012628, -0.45727816, -0.46611571, 0.58052632, -0.24545585, 0.17155244, -0.03935563, 0.64700679, -0.25349086, 0.50070863, 0.45577458, -0.18114147, -0.33747935, 0.64620104, 0.11937550, -0.42006188, -0.45651768;
             0.69254672, -0.03493724, -0.00395996, -0.68322532, -0.27995540, -0.15366832, 0.04854676, 0.51172712, 0.04621797, 0.32861276, 0.04255746, 0.31683190, 0.10598819, -0.45479039, 0.44127824, 0.14584373, -0.46416913, 0.30544929, 0.16581404, -0.61713957, -1.61487831;
             -0.37896436, -0.35747894, -0.48692492, 0.18845077, -0.22663546, 0.15233856, -0.14080742, 0.37993960, 0.28573684, 0.56861321, -0.05663265, -0.52916281, 0.47184619, -0.00726004, -0.49349980, -0.09398657, -0.29787356, -0.39802906, 0.39784611, -0.37318890, -0.70985797;
             0.11630938, 0.25990278, 0.10610731, 0.90468023, -0.72339551, 0.14433450, -0.31667475, 0.19526978, 0.55498026, -0.34111571, -0.11948917, -0.40208843, 0.56907927, -0.19159823, -0.42775992, 0.13477571, -0.11177840, -0.51423039, -0.22485579, -0.06918248, 1.88889082;
             0.19390609, -0.25488125, 0.07153860, -0.55585966, 0.59703165, -0.49822444, -0.04803093, -0.19878294, 0.28714469, -0.45729577, 0.60657821, -0.21444882, 0.18653880, 0.48662769, 0.33005798, -0.07695586, 0.35535236, -0.47732692, -0.08873215, 0.04975119, -0.39763133;
             -0.04516986, -0.48819359, -0.33723572, -0.24124217, 0.44899825, 0.53287724, 0.04801005, 0.04971575, 0.04167457, -0.12421738, 0.48340592, -0.39999628, 0.53383578, 0.27079469, 0.37870918, -0.40744107, -0.07044877, -0.61047126, -0.57235257, 0.30244167, -0.61684750;
             -0.59783935, 0.02848846, -0.30815363, 0.25248297, -0.22056644, 0.61923862, 0.14485995, 0.18447162, -0.32886931, -0.48486822, 0.59139730, -0.25057115, -0.05962549, 0.18348375, -0.28363307, 0.50444079, 0.01209362, -0.33372641, -0.22431492, -0.07320468, -1.21971191],
    HiddenLayer_biasTable = [1.80956654;
             -1.91757424;
             -1.68253726;
             -1.32767956;
             1.25668515;
             -1.15112154;
             -1.61468050;
             0.88897360;
             -0.79853783;
             -0.16027361;
             0.31153760;
             0.29674425;
             -0.13248630;
             -0.45804332;
             0.15665106;
             -0.61836370;
             -0.18421531;
             -1.11280725;
             -1.06571350;
             -1.08889114],
    HiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.TanSig,
    OutputLayer_numNeurons= 1,
    OutputLayer_numInputs= 20,
    OutputLayer_weightTable = [0.17641085, 1.28648182, -0.67012083, 0.70819178, -1.13515224, -1.10506987, 3.50731449, -1.74998461, -0.60570173, 0.59943814, -0.19827507, -1.33285771, 1.22075847, -1.03357541, -1.76913066, 0.12703884, 2.15938712, 0.10077467, 0.36574724, 0.69281302],
    OutputLayer_biasTable = [-0.38211067],
    OutputLayer_NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.LogSig,
        samplePeriod=0.01)
          annotation (extent=[-18,62; 14,90]);

      Modelica.Blocks.Sources.CombiTimeTable out(
        tableOnFile=true,
        tableName="out",
        fileName="testData_RecurrentNN.mat")
                              annotation (extent=[-92,-36; -72,-16]);
      Utilities.SamplerMIMO samplerMIMO(samplePeriod=0.01)
        annotation (extent=[-62,66; -42,86]);
      Utilities.SamplerMIMO samplerMIMO_out(samplePeriod=0.01)
        annotation (extent=[-48,-36; -28,-16]);
      Modelica.Blocks.Sources.CombiTimeTable in2(
        tableOnFile=true,
        tableName="in",
        fileName="testData_RecurrentNN.mat")
                              annotation (extent=[-92,14; -72,34]);
      Networks.NeuralNetwork_RecurrentOneLayer neuralNetwork_RecurrentOneLayer1
        (
        HiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.
            ActivationFunction.TanSig,
        OutputLayer_numNeurons=1,
        HiddenLayer_numInputs=1,
        samplePeriod=0.01,
        HiddenLayer_weightTable=DataFiles.readMATmatrix("testData_RecurrentNN.mat",
            "tanL_weight"),
        HiddenLayer_biasTable=DataFiles.readMATmatrix("testData_RecurrentNN.mat",
            "tanL_bias"),
        OutputLayer_weightTable=DataFiles.readMATmatrix("testData_RecurrentNN.mat",
            "logL_weight"),
        OutputLayer_biasTable=DataFiles.readMATmatrix("testData_RecurrentNN.mat",
            "logL_bias"),
        HiddenLayer_numNeurons=20,
        OutputLayer_numInputs=20,
        OutputLayer_NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.
             LogSig)                   annotation (extent=[-18,10; 14,38]);
      Utilities.SamplerMIMO samplerMIMO2(
                                        samplePeriod=0.01)
        annotation (extent=[-60,14; -40,34]);
    equation
      connect(in1.y, samplerMIMO.u) annotation (points=[-71,76; -64,76], style(
            color=74, rgbcolor={0,0,127}));
      connect(samplerMIMO.y, neuralNetwork_RecurrentOneLayer.u)
                                                               annotation (
          points=[-41,76; -21.2,76], style(color=74, rgbcolor={0,0,127}));
      connect(out.y, samplerMIMO_out.u)
                                     annotation (points=[-71,-26; -50,-26],
          style(color=74, rgbcolor={0,0,127}));
      connect(in2.y, samplerMIMO2.u)
                                    annotation (points=[-71,24; -62,24], style(
            color=74, rgbcolor={0,0,127}));
      connect(samplerMIMO2.y, neuralNetwork_RecurrentOneLayer1.u)
                                                               annotation (
          points=[-39,24; -21.2,24], style(color=74, rgbcolor={0,0,127}));
    end RecurrentNeuralNetwork;

    model RadialBasisNeuralNetwork "An example of a RadialBasis NeuralNetwork"
      annotation (Diagram(
          Rectangle(extent=[-98,98; 96,32], style(color=3, rgbcolor={0,0,255})),
          Text(
            extent=[-84,48; -52,32],
            style(color=3, rgbcolor={0,0,255}),
            string="extractData.m"),
          Rectangle(extent=[-98,26; 96,-40], style(color=3, rgbcolor={0,0,255})),
          Text(
            extent=[-84,-30; -58,-36],
            style(color=3, rgbcolor={0,0,255}),
            string="DataFiles")),
          Documentation(info="<html>
<p>
This model show an example of a RadialBasis NeuralNetwork created with the <i>Neural Network library</i>.
A Radial Basis Neural Network is composed by a layer of neurons which adopts the RadBas activation function and a layer of newrons which perform a purelin transformation.
</p>

<p>
This example is made using a MatLab Radial Basis demo: demorb1. Calling the demo file a new Radial Basis Network (net) is created and trained. The detail about the network created can be see exploring the network object created by matlab. Which is important for our example are these informations:
<UL>
<LI> The first layer of the network, the Radial Basis one, expect 1 input (size of net.IW{1});
<LI> The first layer of the network, the Radial Basis one, is made by 6 newrons (size of net.IW{1});
<LI> The weight table of the first layer is net.IW{1};
<LI> The bias table of the first layer is net.b{1};
<LI> The second layer of the network, the PureLin one, expect 6 inputs (size of net.LW{2});
<LI> The second layer of the network, the PureLin one, is made by 1 newron (size of net.LW{2});
<LI> The weight table of the second layer is net.LW{2};
<LI> The bias table of the second layer is net.b{2};</UL>
</p>

<p>
In this example the parameters of the network are specified using two different way:
<UL>
<LI> using the extractData.m MatLab script, located in Utilities folder;
<LI> using the DataFiles Dymola library.
</UL>
</p>

<p>
The model is simulated on the same data used by demorb1 to test the network created; this can be done creating the file testData_RadialBasisNN using this commands:
<UL>
<LI> call in matlab demorb1
<LI> time=(1:length(X'))'*0.01;
<LI> in=[time,X'];
<LI> out=[time,Y'];
<LI> radL_weight=net.IW{1};
<LI> radL_bias=net.b{1};
<LI> linL_weight=net.LW{2,1};
<LI> linL_bias=net.b{2};
<LI> typing this command in matlab: save testData_RadialBasisNN -V4 in out radL_weight radL_bias linL_weight linL_bias
</UL>
</p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"),
        DymolaStoredErrors);

      Modelica.Blocks.Sources.CombiTimeTable in1(
        tableOnFile=true,
        tableName="in",
        fileName="testData_RadialBasisNN.mat")
                              annotation (extent=[-92,72; -72,92]);
      BaseClasses.NeuralNetworkLayer RadialBasisLayer(
        NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.RadBas,
        numNeurons=6,
        numInputs=1,
        weightTable=[-1.00000000; -0.90000000; -0.80000000; -0.70000000;
            1.00000000; -0.60000000],
        biasTable=[0.83255461; 0.83255461; 0.83255461; 0.83255461; 0.83255461;
            0.83255461])
        annotation (extent=[-52,72; -32,92]);
      BaseClasses.NeuralNetworkLayer PureLinLayer(
        NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.PureLin,
        numNeurons=1,
        numInputs=6,
        weightTable=[99889.17692347,-395573.90463284,593914.14748513,-401122.21457943,
            -123.16770312,102861.39773350],
        biasTable=[85.41689018])
                      annotation (extent=[-12,72; 8,92]);
      Modelica.Blocks.Sources.CombiTimeTable out(
        tableOnFile=true,
        tableName="out",
        fileName="testData_RadialBasisNN.mat")
                              annotation (extent=[-94,-80; -74,-60]);
      Networks.NeuralNetwork_RadialBasis neuralNetwork_RadialBasis(
    HiddenLayer_numNeurons= 6,
    HiddenLayer_numInputs= 1,
    HiddenLayer_weightTable = [-1.00000000;
      -0.90000000;
      -0.80000000;
      -0.70000000;
      1.00000000;
      -0.60000000],
    HiddenLayer_biasTable = [0.83255461;
      0.83255461;
      0.83255461;
      0.83255461;
      0.83255461;
      0.83255461],
    OutputLayer_numNeurons= 1,
    OutputLayer_numInputs= 6,
    OutputLayer_weightTable = [99889.17692347, -395573.90463284, 593914.14748513, -401122.21457943, -123.16770312, 102861.39773350],
    OutputLayer_biasTable = [85.41689018]) annotation (extent=[24,38; 84,82]);
      Modelica.Blocks.Sources.CombiTimeTable in2(
        tableOnFile=true,
        tableName="in",
        fileName="testData_RadialBasisNN.mat")
                              annotation (extent=[-92,0; -72,20]);
      BaseClasses.NeuralNetworkLayer RadialBasisLayer1(
        NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.RadBas,
        numNeurons=6,
        numInputs=1,
        weightTable=DataFiles.readMATmatrix("testData_RadialBasisNN.mat",
            "radL_weight"),
        biasTable=DataFiles.readMATmatrix("testData_RadialBasisNN.mat",
            "radL_bias"))
        annotation (extent=[-52,0; -32,20]);
      BaseClasses.NeuralNetworkLayer PureLinLayer1(
        NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.PureLin,
        numNeurons=1,
        numInputs=6,
        weightTable=DataFiles.readMATmatrix("testData_RadialBasisNN.mat",
            "linL_weight"),
        biasTable=DataFiles.readMATmatrix("testData_RadialBasisNN.mat",
            "linL_bias"))
                      annotation (extent=[-12,0; 8,20]);
      Networks.NeuralNetwork_RadialBasis neuralNetwork_RadialBasis1(
    HiddenLayer_numNeurons= 6,
    HiddenLayer_numInputs= 1,
    OutputLayer_numNeurons= 1,
    OutputLayer_numInputs= 6,
        HiddenLayer_weightTable=DataFiles.readMATmatrix(
            "testData_RadialBasisNN.mat", "radL_weight"),
        HiddenLayer_biasTable=DataFiles.readMATmatrix(
            "testData_RadialBasisNN.mat", "radL_bias"),
        OutputLayer_weightTable=DataFiles.readMATmatrix(
            "testData_RadialBasisNN.mat", "linL_weight"),
        OutputLayer_biasTable=DataFiles.readMATmatrix(
            "testData_RadialBasisNN.mat", "linL_bias"))
        annotation (extent=[24,-34; 84,10]);
    equation
      connect(RadialBasisLayer.y, PureLinLayer.u)
        annotation (points=[-31,82; -14,82],
                                          style(color=74, rgbcolor={0,0,127}));
      connect(in1.y, RadialBasisLayer.u) annotation (points=[-71,82; -54,82], style(
            color=74, rgbcolor={0,0,127}));
      connect(in1.y, neuralNetwork_RadialBasis.u) annotation (points=[-71,82;
            -66,82; -66,60; 18,60], style(color=74, rgbcolor={0,0,127}));
      connect(RadialBasisLayer1.y, PureLinLayer1.u)
        annotation (points=[-31,10; -14,10],
                                          style(color=74, rgbcolor={0,0,127}));
      connect(in2.y, RadialBasisLayer1.u)
                                         annotation (points=[-71,10; -54,10], style(
            color=74, rgbcolor={0,0,127}));
      connect(in2.y, neuralNetwork_RadialBasis1.u) annotation (points=[-71,10;
            -66,10; -66,-12; 18,-12], style(color=74, rgbcolor={0,0,127}));
    end RadialBasisNeuralNetwork;
    annotation (Documentation(info="<html>

In this package there are some examples of using the NeuralNetwork library.

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"));
    model FeedForwardNeuralNetwork
      "An example of a Feed-forward NeuralNetwork: in particultar, it was implemented using Networks elements"

      Modelica.Blocks.Discrete.UnitDelay unitDelay(samplePeriod=0.01)
        annotation (extent=[-54,36; -48,42]);
      Modelica.Blocks.Routing.Multiplex3 multiplex3_1
        annotation (extent=[-26,20; -20,56]);
      annotation (Diagram, Icon,
        Documentation(info="<html>
<p>
This model show an example of a Feed-forward NeuralNetwork created with the <i>Neural Network library</i>. This is composed by two layer: the hidden one uses a TanSig activation function and the output one uses a PureLin activation function.
</p>


<p>
This example is made using the following MatLab commands:
<UL>
<LI> clear
<LI> t=0:0.01:10;
<LI> x=sin(2*pi*t);
<LI> y=cos(5*pi*t);
<LI> for k=3:length(t) f(k)=(3*x(k)*x(k-1)*x(k-2))+(y(k)*y(k-2)); end
<LI> var_X = [[min(x) max(x)];[min(x) max(x)];[min(x) max(x)]];
<LI> var_Y = [[min(y) max(y)];[min(y) max(y)];[min(y) max(y)]];
<LI> in_X = [x ; [ 0 x(1:end-1)] ; [ 0 0 x(1:end-2)]];
<LI> in_Y = [y ; [ 0 y(1:end-1)] ; [ 0 0 y(1:end-2)]];
<LI> net = newff([var_X ; var_Y],[4 1],{'tansig','purelin'});
<LI> net.trainFcn = 'trainlm';
<LI> net.trainParam.epochs = 100;
<LI> [net,tr] = train(net , [in_X ; in_Y] , f);
<LI> f_SIM = sim(net,[in_X;in_Y]);
<LI> figure;
<LI> hold;
<LI> plot(t,f);
<LI> plot(t,f_SIM,'r');
<LI> IN_x=[t' , x'];
<LI> IN_y=[t' , y'];
<LI> OUT_f=[t' , f_SIM'];
<LI> save testData_FeedForwardNN.mat -V4 IN_x IN_y OUT_f
<LI> EXTRACTDATA('LW.txt','OutputLayer',net.LW{2,1},net.b{2},'lin')
<LI> EXTRACTDATA('IW.txt','HiddenLayer',net.IW{1},net.b{1},'tan')
</UL>

These commands create and train a new FeedForward Neural Network (net), which takes in input x and y and them delayed. The detail about the network created can be see exploring the network object created by matlab. Which is important for our example are these informations:
<UL>
<LI> The first layer of the network, the TanSig one, expect 6 input (size of net.IW{1});
<LI> The first layer of the network, the TanSig one, is made by 4 newrons (size of net.IW{1});
<LI> The weight table of the first layer is net.IW{1};
<LI> The bias table of the first layer is net.b{1};
<LI> The second layer of the network, the LogSig one, expect 4 inputs (size of net.LW{2});
<LI> The second layer of the network, the LogSig one, is made by 1 newron (size of net.LW{2});
<LI> The weight table of the second layer is net.LW{2};
<LI> The bias table of the second layer is net.b{2};</UL>
</p>

<p>
In this example the parameters of the network are specified using the extractData.m MatLab script, located in Utilities folder.
</p>

<p>
The model is simulated on the same data used in MatLab so it is possible to test if the model gives the expected output.
<p>

<p><b>Release Notes:</b></p>
<ul>
<li><i>April 28, 2006</i>
       by <a href=\"http://home.dei.polimi.it/codeca\">Fabio Codecà</a>:<br></li>
</ul>
</HTML>"));
      Modelica.Blocks.Discrete.UnitDelay unitDelay6(samplePeriod=0.01)
        annotation (extent=[-54,-24; -48,-18]);
      Modelica.Blocks.Routing.Multiplex3 multiplex3_2
        annotation (extent=[-26,-38; -20,-2]);
      Modelica.Blocks.Discrete.UnitDelay unitDelay1(samplePeriod=0.01)
        annotation (extent=[-42,24; -36,30]);
      Modelica.Blocks.Discrete.UnitDelay unitDelay3(samplePeriod=0.01)
        annotation (extent=[-44,-34; -38,-28]);
      Modelica.Blocks.Discrete.Sampler sampler1(samplePeriod=0.01)
        annotation (extent=[-72,-10; -64,-2]);
      Modelica.Blocks.Discrete.Sampler sampler3(samplePeriod=0.01)
        annotation (extent=[-70,48; -62,56]);

      Networks.NeuralNetwork_TwoLayer neuralNetwork_TwoLayer1(
       HiddenLayer_numNeurons= 4,
    HiddenLayer_numInputs= 6,
    HiddenLayer_weightTable = [0.42245576, -0.25417188, 0.39975683, -0.07080193, 0.04146637, -0.17979616;
      0.41799358, -0.86810195, 0.41350926, -0.25591464, 0.27097464, -0.50095646;
      0.48053988, -0.50621676, 0.48146247, 0.04731258, -0.16416807, -0.02218562;
      1.01725958, -0.86088932, 1.00258563, 0.24379596, -0.48057624, 0.24876107],
    HiddenLayer_biasTable = [-0.86790396;
      0.69253240;
      0.59726118;
      0.00460979],
    HiddenLayer_NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.TanSig,OutputLayer_numNeurons= 1,
    OutputLayer_numInputs= 4,
    OutputLayer_weightTable = [13.86326791, -11.79653122, 20.90177696, -10.02653340],
    OutputLayer_biasTable = [5.59944933],
    OutputLayer_NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.PureLin)
                                        annotation (extent=[40,0; 60,20]);
    public
      Modelica.Blocks.Routing.Multiplex2 multiplex2_1(n1=3, n2=3)
        annotation (extent=[4,0; 24,20]);
      Modelica.Blocks.Sources.CombiTimeTable IN_x(
        tableOnFile=true,
        fileName="testData_FeedForwardNN.mat",
        tableName="IN_x")     annotation (extent=[-100,42; -80,62]);
      Modelica.Blocks.Sources.CombiTimeTable IN_y(
        tableOnFile=true,
        fileName="testData_FeedForwardNN.mat",
        tableName="IN_y")     annotation (extent=[-100,-16; -80,4]);
      Modelica.Blocks.Sources.CombiTimeTable OUT_f(
        tableOnFile=true,
        fileName="testData_FeedForwardNN.mat",
        tableName="OUT_f")    annotation (extent=[-100,-68; -80,-48]);
      Modelica.Blocks.Discrete.Sampler sampler_OUT(samplePeriod=0.01)
        annotation (extent=[-72,-62; -64,-54]);
    equation
      connect(unitDelay6.y, multiplex3_2.u2[1]) annotation (points=[-47.7,-21;
            -38.85,-21; -38.85,-20; -26.6,-20],
                                          style(color=74, rgbcolor={0,0,127}));
      connect(unitDelay.y, multiplex3_1.u2[1]) annotation (points=[-47.7,39;
            -35.85,39; -35.85,38; -26.6,38], style(color=74, rgbcolor={0,0,127}));
      connect(unitDelay.y, unitDelay1.u) annotation (points=[-47.7,39; -46,34;
            -46,28; -42,28; -42.6,27], style(color=74, rgbcolor={0,0,127}));
      connect(unitDelay1.y, multiplex3_1.u3[1]) annotation (points=[-35.7,27;
            -31.85,27; -31.85,25.4; -26.6,25.4], style(color=74, rgbcolor={0,0,
              127}));
      connect(unitDelay6.y, unitDelay3.u) annotation (points=[-47.7,-21; -47.7,
            -26.5; -44.6,-26.5; -44.6,-31],
                                       style(color=74, rgbcolor={0,0,127}));
      connect(unitDelay3.y, multiplex3_2.u3[1]) annotation (points=[-37.7,-31;
            -33.85,-31; -33.85,-32.6; -26.6,-32.6],
                                                 style(color=74, rgbcolor={0,0,
              127}));
      connect(sampler1.y, multiplex3_2.u1[1]) annotation (points=[-63.6,-6;
            -35.1,-6; -35.1,-7.4; -26.6,-7.4], style(color=74, rgbcolor={0,0,
              127}));
      connect(sampler1.y, unitDelay6.u) annotation (points=[-63.6,-6; -58,-6;
            -58,-21; -54.6,-21],
                             style(color=74, rgbcolor={0,0,127}));
      connect(sampler3.y, multiplex3_1.u1[1]) annotation (points=[-61.6,52;
            -44.1,52; -44.1,50.6; -26.6,50.6], style(color=74, rgbcolor={0,0,
              127}));
      connect(sampler3.y, unitDelay.u) annotation (points=[-61.6,52; -58,52;
            -58,39; -54.6,39], style(color=74, rgbcolor={0,0,127}));
      connect(multiplex3_1.y, multiplex2_1.u1) annotation (points=[-19.7,38;
            -14,38; -14,16; 2,16], style(color=74, rgbcolor={0,0,127}));
      connect(multiplex3_2.y, multiplex2_1.u2) annotation (points=[-19.7,-20;
            -14,-20; -14,4; 2,4], style(color=74, rgbcolor={0,0,127}));
      connect(multiplex2_1.y, neuralNetwork_TwoLayer1.u) annotation (points=[25,10;
            38,10],     style(color=74, rgbcolor={0,0,127}));
      connect(IN_x.y[1], sampler3.u) annotation (points=[-79,52; -70.8,52],
          style(color=74, rgbcolor={0,0,127}));
      connect(IN_y.y[1], sampler1.u) annotation (points=[-79,-6; -72.8,-6],
          style(color=74, rgbcolor={0,0,127}));
      connect(OUT_f.y[1], sampler_OUT.u)
                                      annotation (points=[-79,-58; -72.8,-58],
          style(color=74, rgbcolor={0,0,127}));
    end FeedForwardNeuralNetwork;
  end Examples;
end NeuralNetwork;
