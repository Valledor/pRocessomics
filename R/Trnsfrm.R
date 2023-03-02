#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 4.11.2019
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


Trnsfrm <- function(matriz, transformation="z", initialrow=1, initialcolumn=3,treatment1col=1, 
                    treatment2col=2, treatment=1,datasetname=NULL){
  cat(paste("\nProcessing ",datasetname," dataset",sep=""))
  datos<-matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)]
  switch(transformation,
         "z"=datostransformados <- scale(datos,center=F),
         "z-centered"=datostransformados <- scale(datos),
         "Log10"={ceros<-datos==0
         datostransformados <- log(datos,10)
         datostransformados[ceros==T]<-0},
         "Log10plus1"=datostransformados<-log(datos+1,10),
         "box-cox-lambda"={
           if(treatment==1) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment1col])
           if(treatment==2) vectortratamientos <- (matriz[initialrow:nrow(matriz),treatment2col])
           if(treatment==3) vectortratamientos <- (paste(matriz[initialrow:nrow(matriz),treatment1col],matriz[initialrow:nrow(matriz),treatment2col], sep="&"))
           vectortratamientos<-factor(vectortratamientos, levels = unique(vectortratamientos))
           datostransformados<-apply(as.data.frame(datos),2,function(x) boxcoxvariable(x,vectortratamientos))
         },
         "sqrt"=datostransformados <- datos^(1/2),
         "cubic"=datostransformados <- datos^(1/3),
         "arcsine"=datostransformados <- asin(datos/100)*2/pi,
         "none"=datostransformados<-datos,
         return(cat("\nTransformation not recognized. Please use 'z','z-centered','Log10', 'log10plus1', 'sqrt', 'cubic', 'arcsine', 'box-cox-lambda' or 'none'.")))
  matriz[initialrow:nrow(matriz),initialcolumn:ncol(matriz)] <-datostransformados
  return(matriz)
}
