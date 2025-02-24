#
# Copyright 2019 Luis Valledor / Laura Lamelas
# Last revision 04.11.2019
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
feature_selection_list_helper <- function(tabladatos, initialrow, initialcolumn, treatment1col, treatment2col, treatment, method, threshold, parametric,stat,datasetname=NULL,abovebelow){
  if(method=="IQR") tabladatos<-IQRfilter(tabladatos, initialrow = initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, treatment=treatment, threshold=threshold,abovebelow=abovebelow)
  if(method=="Stat") tabladatos<-Statfilter(tabladatos, initialrow = initialrow, initialcolumn=initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, treatment = treatment, threshold=threshold,parametric=parametric, stat=stat,abovebelow=abovebelow)
  if(method=="VarCoef") tabladatos<-VarCoefSelect(tabladatos, initialrow = initialrow, initialcolumn = initialcolumn, treatment1col=treatment1col, treatment2col=treatment2col, treatment = treatment, threshold=threshold, datasetname = datasetname,abovebelow=abovebelow)
  if(method=="none") {
    tabladatos<-tabladatos
    texts_wizard("\nDoing nothing!")
    }
  return(tabladatos)
}
