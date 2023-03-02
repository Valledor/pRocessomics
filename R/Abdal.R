#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 04.02.2020
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


AbdBal <-
  function(matriz,
           norm = "AvgIntensity",
           initialrow = 1,
           initialcolumn = 3,
           treatment1col = 1,
           treatment2col = 2,
           treatment = 1,
           datasetname = NULL) {
    cat(paste("\nProcessing ", datasetname, " dataset", sep = ""))
    if (norm == "none") {
      return(matriz)
    }
    #### Sample ####
    else if (norm == "Sample") {
      datos <- matriz[initialrow:nrow(matriz), initialcolumn:ncol(matriz)]
      matriznormalizada <-
        t(apply(datos, 1, function(x)
          x / sum(x, na.rm = T)))
      matriz[initialrow:nrow(matriz), initialcolumn:ncol(matriz)] <-
        matriznormalizada
      return(matriz)
    }
    #### AvgIntensity ####
    else if (norm == "AvgIntensity") {
      datos <- matriz[initialrow:nrow(matriz), initialcolumn:ncol(matriz)]
      abundanciamedia <- mean(rowSums(datos, na.rm = T), na.rm = T)
      matriznormalizada <-
        t(apply(datos, 1, function(x)
          x * abundanciamedia / sum(x, na.rm = T)))
      matriz[initialrow:nrow(matriz), initialcolumn:ncol(matriz)] <-
        matriznormalizada
      return(matriz)
    }
    #### TreatAvgIntensity ####
    else if (norm == "TreatAvgIntensity") {
      if (treatment == 1)
        vectortratamientos <-
          (matriz[initialrow:nrow(matriz), treatment1col])
      if (treatment == 2)
        vectortratamientos <-
          (matriz[initialrow:nrow(matriz), treatment2col])
      if (treatment == 3)
        vectortratamientos <-
          (paste(matriz[initialrow:nrow(matriz), treatment1col], matriz[initialrow:nrow(matriz), treatment2col], sep =
                   "&"))
      vectortratamientos <-
        factor(vectortratamientos, levels = unique(vectortratamientos))
      
      datos <- split(as.data.frame(matriz), vectortratamientos, drop = T)
      processedlist <-
        lapply(datos, function(x)
          AbdBalTreatment(x, initialcolumn))
      matriz2 <- do.call("rbind",  unname(processedlist))
      matriz2 <- matriz2[match(rownames(matriz), rownames(matriz2)),]
      return(matriz2)
    }
  }
