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

#' @importFrom MASS boxcox

boxcoxvariable <-function(vector,tratamiento){
  #suppressMessages(suppressWarnings(require(MASS)))
  ceros <- vector==0 
  vector[vector==0] <- 0.000000000001
  lambda <- seq(-5,5,0.2)
  box <- MASS::boxcox(vector~as.factor(tratamiento),lambda,plotit=F) 
  cox <- data.frame(box$x, box$y) 
  cox2 <- cox[with(cox, order(-cox$box.y)),] 
  lambda <- cox2[1,"box.x"] 
  if(lambda ==0){
    t_vector <-log(vector)
  } else{
    t_vector <- (vector^lambda)#-1/lambda   #(vector^lambda)-1/lambda Forma original de BoxCox
  }
  t_vector[ceros==T] <-0
  return(t_vector)
}
