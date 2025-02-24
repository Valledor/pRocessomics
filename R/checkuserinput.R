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


checkuserinput <-function(variable,nrow=NULL,ncol=NULL){
  if(is.numeric(variable)==F) stop("FATAL ERROR. Please input an integer number corresponding to desired row/column.")
  if(variable%%1!=0) stop("FATAL ERROR.Please select a correct row/column by introducing an integer. Please re-run.")
  if(is.null(nrow)) if(variable > ncol) stop("FATAL ERROR. The column you have selected is not available in your dataset. Please re-run and select a correct column.")
  if(is.null(ncol)) if(variable > nrow) stop("FATAL ERROR. The row you have selected is not available in your dataset. Please re-run and select a correct row.")
}
