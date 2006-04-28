%EXTRACTDATA Create parameters code for Dymola NeuralNetwork Library.
%
%   EXTRACTDATA(FILENAME,LAYERNAME,MATRIXWEIGHT,MATRIXBIAS,ACTFUN,RECURRENT)
%   creates a file named as FILENAME, in which is created the code for a
%   NeuronNetwork layer.
%
%   If the NewronNetworkLayer is used the LAYERNAME has to be equal to
%   'NoLayerName', otherwise as to be equal to the internal layer which
%   parameters are referred.
%
%   If the script is used for a recurrent layer the parameter RECURRENT has
%   to be equal to the number of non-recurrent inputs to the layer. For
%   example, if the layer has 1 non-recurrent input and 5 neurons then the
%   all input to the layer will be 6, but RECURRENT has to be equal to 1.

%   ACTFUN can be equal to:
%       'lin' for PureLin NeuralNetwork Activation Function
%       'tan' for TanSig NeuralNetwork Activation Function
%       'log' for LogSig NeuralNetwork Activation Function
%       'rad' for RadBas NeuralNetwork Activation Function
%       'noFun' is the activation function has not to be specified
%
%   Made by Fabio Codecà (http://www.elet.polimi.it/people/codeca)



function extractData(fileName,layerName,matrixWeight,matrixBias,actFun,recurrent)

if (strcmp('NoLayerName',layerName))
    layerName='';
else
    layerName=strcat(layerName,'_');
end

fid = fopen(fileName,'w');
Y = [matrixWeight];
Z = [matrixBias];
[numRow, numCol] = size(Y)

if nargin == 5
    numInput = numCol;
else
    numInput = recurrent;
end

fprintf(fid,strcat(layerName,'numNeurons = ',num2str(numRow),', \n'));
fprintf(fid,strcat(layerName,'numInputs = ',num2str(numInput),', \n'));

fprintf(fid,strcat(layerName,'weightTable = ['));
% j row
% i column
for j=1:numRow
    if (numCol==1)
        if (numRow == 1)
            fprintf(fid,'%8.8f ], \n',Y(j,1));
        elseif (j == 1)
            fprintf(fid,'%8.8f ; \n',Y(j,1));
        elseif (j < numRow)
            fprintf(fid,'\t %8.8f ; \n',Y(j,1));
        else
            fprintf(fid,'\t %8.8f ],\n',Y(j,1));
        end
    else
        for i=1:numCol
            if (j > 1 && i == 1)
                fprintf(fid,'\t %8.8f ,',Y(j,i));
            elseif (i==numCol)
                if (j < numRow)
                    fprintf(fid,'%8.8f ; \n',Y(j,i));
                else
                    fprintf(fid,'%8.8f ],\n',Y(j,i));
                end
            else
                fprintf(fid,'%8.8f ,',Y(j,i));
            end
        end
    end
end

numCol=1;
fprintf(fid,strcat(layerName,'biasTable = ['));
for j=1:numRow
    if (numCol==1)
        if (numRow == 1)
            fprintf(fid,'%8.8f ]',Z(j,1));
        elseif (j == 1)
            fprintf(fid,'%8.8f ; \n',Z(j,1));
        elseif (j < numRow)
            fprintf(fid,'\t %8.8f ; \n',Z(j,1));
        else
            fprintf(fid,'\t %8.8f ]',Z(j,1));
        end
    else
        for i=1:numCol
            if (j > 1 && i == 1)
                fprintf(fid,'\t %8.8f ,',Z(j,i));
            elseif (i==numCol)
                if (j < numRow)
                    fprintf(fid,'%8.8f ; \n',Z(j,i));
                else
                    fprintf(fid,'%8.8f ]',Z(j,i));
                end
            else
                fprintf(fid,'%8.8f ,',Z(j,i));
            end
        end
    end
end


printFun=1;

if (strcmp(actFun,'lin'))
    fun='PureLin';
elseif (strcmp(actFun,'tan'))
    fun='TanSig';
elseif (strcmp(actFun,'log'))
    fun='LogSig';
elseif (strcmp(actFun,'rad'))
    fun='RadBas';
elseif (strcmp(actFun,'noFun'))
    printFun = 0;
else
    error('Activation Function Name is incorrect');
end

if(printFun)
    fprintf(fid,strcat(', \n',layerName,'NeuronActivationFunction=NeuralNetwork.Types.ActivationFunction.',fun));
end

fclose(fid);
