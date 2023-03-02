##
## Definimos paletas y simbolos para para plots. Si queremos dar opcion al usuario de escoger
## varias paletas, esto se puede sacar fuera de la funcion y pasar las paletas como argumentos
##

#### Default color palettes and vectors ####
paleta_mapman <- c("#59c260", "#55823b", "#bcdf8d", "#dadd42", "#88bc3a", "#79e54a", 
                   "#a6a142", "#c9cfa6", "#cb8238", "#BDB76B", "#FFD700", "#d54837", 
                   "#c07d77", "#ce4867", "#e697a4", "#e53dc8", "#d43f6f", "#9b5b75", 
                   "#cd57a3", "#e391c8", "#e391c8", "#656f91", "#7f53e0", "#e391c8",
                   "#8896dc", "#DCDCDC", "#e391c8", "#a073b9", "#4861ac", "#6bccd5",
                   "#4c7a95", "#4c7a95", "#74b0e0", "#78a7a7", "#808080")
paleta_tratamientos <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231","#c217b8", 
                         "#38d1d1", "#c24f98")
paleta_niveles <- c("#de6ab7","#ec4e20","#3a7ca5","#1e8a7d","#d4aa13","#9e470d",
                    "#691e0c","#4f2828")
paleta_continuo <- c("#5E4FA2", "#5954A4", "#555AA7", "#5060AA", "#4C66AC", "#476BAF", 
                     "#4371B2", "#3E77B5", "#3A7DB7", "#3682BA", "#3288BC", "#378EBA", 
                     "#3D94B7", "#429AB5", "#47A0B3", "#4CA5B0", "#52ABAE", "#57B1AB", 
                     "#5CB7A9", "#61BDA6", "#67C2A4", "#6EC5A4", "#75C8A4", "#7CCAA4", 
                     "#83CDA4", "#8AD0A4", "#91D2A4", "#98D5A4", "#9FD8A4", "#A6DBA4", 
                     "#ACDDA3", "#B2E0A2", "#B8E2A1", "#BEE5A0", "#C4E79E", "#CAE99D", 
                     "#D0EC9C", "#D6EE9B", "#DCF199", "#E2F398", "#E7F599", "#E9F69D", 
                     "#ECF7A1", "#EEF8A5", "#F1F9A9", "#F3FAAD", "#F6FBB1", "#F8FCB5", 
                     "#FBFDB9", "#FDFEBD", "#FEFDBC", "#FEFAB7", "#FEF7B1", "#FEF4AC", 
                     "#FEF0A7", "#FEEDA2", "#FEEA9C", "#FEE797", "#FEE492", "#FEE18D", 
                     "#FDDC88", "#FDD784", "#FDD27F", "#FDCD7B", "#FDC877", "#FDC373", 
                     "#FDBE6F", "#FDB96A", "#FDB466", "#FDAF62", "#FCA95E", "#FBA25B", 
                     "#FA9C58", "#F99555", "#F88F52", "#F7884F", "#F6824C", "#F67B49", 
                     "#F57446", "#F46E43", "#F16943", "#EE6445", "#EB5F46", "#E85A47", 
                     "#E45648", "#E1514A", "#DE4C4B", "#DB474C", "#D8434D", "#D53E4E", 
                     "#CF384D", "#CA324C", "#C42C4B", "#BF2549", "#B91F48", "#B41947", 
                     "#AE1345", "#A90D44", "#A30743", "#9E0142")

# Symbols
simbolos <- c("circle","cross","square","star","diamond","triangle-up")

# Default fontsizes
fontsizes <- c(14,10,16,12)

# Default vectors
truefalse <- c("Yes","No")

#### Function for introducing used-defined palettes ####

#' @name set_custom_palettes
#' @title Set custom palette for pRocessomics plots
#' @description This function allows to create a new palette to be displayed in pRocessomics plots. Including symbols, fontsizes and colors
#' @usage set_custom_palettes(palette_type,newcustomvalues)
#' @param palette_type
#' \itemize{
#' \item "mapman" For Mapman (or Custom annotation) related plots.
#' \item "levels" For different omic levels 
#' \item "treatment" For different treatments
#' \item "continuous" For continuous scales
#' \item "symbols" 
#' \item "fontsizes"
#' }
#' @param newcustomvalues character vector especifying desired values.
#' @return Override the default values to the user specified 
#' @author Luis Valledor and Laura Lamelas
#' 
#' @export
#' @importFrom graphics rect text par plot

#' 
set_custom_palettes <- function(palette_type,newcustomvalues){
  
  if(palette_type %in% c("mapman","levels","treatment","continuous","symbols","fontsizes")==F) stop("FATAL Error: select adequate palette name to replace. Available palettes are \"mapman\",\"levels\",\"treatment\",\"continuous\",\"symbols\",\"fontsizes\"")
  
  if(palette_type =="mapman" & length(newcustomvalues)!=35) stop("FATAL Error: mapman palette should have exactly 35 elements")
  if(palette_type =="continuous" & length(newcustomvalues)!=100) stop("FATAL Error: continuous palette should have exactly 100 elements")
  if(palette_type =="fontsizes" & length(fontsizes)!=4) stop("FATAL Error: font size vector palette should have exactly 4 elements.")
  
  if(palette_type =="mapman") {
    new_mapman_palette <- newcustomvalues
    class(new_mapman_palette) <- c("color.pRo", "mapman")
    return(new_mapman_palette)}
  
  if(palette_type =="levels") {
    new_level_palette <- newcustomvalues
    class(new_level_palette) <- c("color.pRo", "levels")
    return(new_level_palette)}
  
  if(palette_type =="treatment") {
    mew_treatment_palette <- newcustomvalues
    class(mew_treatment_palette) <- c("color.pRo", "treatment")
    return(mew_treatment_palette)}
  
  if(palette_type =="continuous") {
    new_continous_palette <- newcustomvalues
    class(new_continous_palette) <- c("color.pRo", "continous")
    return(new_continous_palette)}
  
  if(palette_type =="symbols") {
    new_symbols <- newcustomvalues
    class(new_symbols) <- c("color.pRo", "symbols")
    return(new_symbols)}
  
  if(palette_type =="fontsizes") {
    new_fontsizes <- newcustomvalues
    class(new_fontsizes) <- c("color.pRo", "fontsizes")
    return(new_fontsizes)}
  
  
}


#' @name select_palette
#' @title Select and store one of the pre-defined palettes for pRocessomics plots, pRocessomics functions will use the selected palettes to display the plots automatically if any pre-selected palette is in the global environment
#' @description This function allows to select de desired palette (among different options) to be displayed in pRocessomics plots. 
#' @usage select_palette()
#' @return charater vector, class "color.pRo" to set user defined colors/symbols/fontsizes to pRocessomics generated plots
#' @author Luis Valledor and Laura Lamelas
#' 
#' @export
#' @importFrom graphics rect text par plot
#' @importFrom utils select.list
#' @importFrom grDevices dev.cur

# Function for selecting a pre-defined palette
select_palette <- function(){
  
  headings_wizards("pRocessomcis - Palette Selection")
  texts_wizard("\n\nImage shows pRocessomics color space. In the upper rows you can see the palettes currently used for levels (Current-Lvls). Pre-designed palettes are shown below. You can also employ your own colors and symbols palettes and import them to pRocessomics with set_custom_palettes(). Please read 'https://github.com/Valledor/pRocessomics/wiki/Pre-processing-your-data' before continuing.")
  texts_wizard("\nWe recommend employing RColorBrewer or Viridis packages for getting beautiful palettes. Part of pRocessomics palettes have been generated with these tools.\n")
  texts_wizard("Press [enter] to continue...")
  invisible(readline())
  # Reset environment
  while (grDevices::dev.cur()>1) {
  dev.off()}
  graphics::par(mfrow=c(1,1),mar=c(0, 4.5, 1.5, 0))
  
  #### Define Palettes ####  
  Colorful <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231","#c217b8", "#38d1d1")
  Suede <- c("#A7226E","#EC2049","#F26B38", "#F7DB4F", "#2F9599", "#547980", "#594F4F")
  Strong <- c("#e11f25", "#3b7eb6", "#51ae4f", "#9851a1", "#fc8023", "#f5f546", "#693417")
  Rainbow <- c("#e6235b","#f08232","#fae650","#a0c84b","#4bbea0","#4182be","#644ba0")
  Calid <- c("#de6ab7","#ec4e20","#3a7ca5","#1e8a7d","#d4aa13","#9e470d","#691e0c")
  Intense <- c("red","blue","green","purple","orange","yellow","brown")
  
  #### Plot colors for omic levels and treatments ####
  listapaletas <- list(Intense, Rainbow, Calid, Strong, Suede, Colorful,paleta_niveles,paleta_tratamientos)
  names(listapaletas) <- c("Intense", "Rainbow", "Calid", "Strong", "Suede", "Colorful", "Current-Lvls", "Current-Trmnt")
  
  n <-sapply(listapaletas,function(x) length(x))
  nr <- length(listapaletas)
  nc <- max(n)
  
  ylim <- c(0,nr)
  oldpar <- graphics::par(mgp=c(2,0.25,0))
  on.exit(graphics::par(oldpar))
  graphics::plot(1,1,xlim=c(0,nc),ylim=ylim,type="n", axes=FALSE, bty="n",
       xlab="",ylab="",main = "pRocessomics color space")
  for(i in 1:nr)
  {nj <- n[i]
  #if (colorlist[i]=="") next
  shadi <- listapaletas[[i]]
  graphics::rect(xleft=0:(nj-1), ybottom=i-1, xright=1:nj, ytop=i-0.2, col=shadi,
       border="light grey")
  }
  graphics::text(rep(-0.2,nr),(1:nr)-0.6, labels=names(listapaletas), xpd=TRUE, adj=1)
  
  # Selecting palette for treatments
  palettes <- c("Current", "Colorful", "Suede", "Strong", "Calid", "Rainbow","Intense")
  choicetreat <- utils::select.list(choices = palettes,multiple = F,preselect = "Current",title = "\nPlease select a number of palette for assigning colors to the different treatments and press [enter].")
  if(choicetreat=="Colorful") new_paleta_tratamientos <- Colorful
  if(choicetreat=="Suede") new_paleta_tratamientos <- Suede
  if(choicetreat=="Strong") new_paleta_tratamientos <- Strong
  if(choicetreat=="Calid") new_paleta_tratamientos <- Calid
  if(choicetreat=="Rainbow") new_paleta_tratamientos <- Rainbow
  if(choicetreat=="Intense") new_paleta_tratamientos <- Intense
  revlvlcols <- utils::select.list(choices = truefalse,multiple = F, title = "\nDo you want to reverse the order of the colors?")
  if(revlvlcols=="Yes") new_paleta_tratamientos<-rev(new_paleta_tratamientos)
  new_treatment_palette <- new_paleta_tratamientos
  class(new_treatment_palette) <- c("color.pRo", "treatment")
  new <- list(new_treatment_palette)
  names(new) <- c("new_treatment_palette")
  list2env(new, envir = .GlobalEnv)
  
  # Selecting palette for omic levels
  choicelvl <- utils::select.list(choices = palettes,multiple = F,preselect = "Current",title = "\nPlease select a number of palette for assigning colors to the different omic levels and press [enter].")
  if(choicelvl=="Colorful") new_paleta_niveles <- Colorful
  if(choicelvl=="Suede") new_paleta_niveles <- Suede
  if(choicelvl=="Strong") new_paleta_niveles <- Strong
  if(choicelvl=="Calid") new_paleta_niveles <- Calid
  if(choicelvl=="Rainbow") new_paleta_niveles <- Rainbow
  if(choicelvl=="Intense") new_paleta_niveles <- Intense
  revtreatcols <- utils::select.list(choices = truefalse,multiple = F, title = "\nDo you want to reverse the order of the colors?")
  if(revtreatcols=="Yes") new_paleta_niveles<-rev(new_paleta_niveles)
  
  new_level_palette <- new_paleta_niveles
  class(new_level_palette) <- c("color.pRo", "levels")
  new <- list(new_level_palette)
  names(new) <- c("new_level_palette")
  list2env(new, envir = .GlobalEnv)
  
  #### Define Continuous paletete ####
  paleta_continuo_1 <- c("#5E4FA2", "#5954A4", "#555AA7", "#5060AA", "#4C66AC", "#476BAF", 
                         "#4371B2", "#3E77B5", "#3A7DB7", "#3682BA", "#3288BC", "#378EBA", 
                         "#3D94B7", "#429AB5", "#47A0B3", "#4CA5B0", "#52ABAE", "#57B1AB", 
                         "#5CB7A9", "#61BDA6", "#67C2A4", "#6EC5A4", "#75C8A4", "#7CCAA4", 
                         "#83CDA4", "#8AD0A4", "#91D2A4", "#98D5A4", "#9FD8A4", "#A6DBA4", 
                         "#ACDDA3", "#B2E0A2", "#B8E2A1", "#BEE5A0", "#C4E79E", "#CAE99D", 
                         "#D0EC9C", "#D6EE9B", "#DCF199", "#E2F398", "#E7F599", "#E9F69D", 
                         "#ECF7A1", "#EEF8A5", "#F1F9A9", "#F3FAAD", "#F6FBB1", "#F8FCB5", 
                         "#FBFDB9", "#FDFEBD", "#FEFDBC", "#FEFAB7", "#FEF7B1", "#FEF4AC", 
                         "#FEF0A7", "#FEEDA2", "#FEEA9C", "#FEE797", "#FEE492", "#FEE18D", 
                         "#FDDC88", "#FDD784", "#FDD27F", "#FDCD7B", "#FDC877", "#FDC373", 
                         "#FDBE6F", "#FDB96A", "#FDB466", "#FDAF62", "#FCA95E", "#FBA25B", 
                         "#FA9C58", "#F99555", "#F88F52", "#F7884F", "#F6824C", "#F67B49", 
                         "#F57446", "#F46E43", "#F16943", "#EE6445", "#EB5F46", "#E85A47", 
                         "#E45648", "#E1514A", "#DE4C4B", "#DB474C", "#D8434D", "#D53E4E", 
                         "#CF384D", "#CA324C", "#C42C4B", "#BF2549", "#B91F48", "#B41947", 
                         "#AE1345", "#A90D44", "#A30743", "#9E0142") 
  
  paleta_continuo_2 <- c("#006837", "#026C39", "#05713C", "#07763E", "#0A7B41", "#0D8043", 
                         "#0F8546", "#128948", "#158E4B", "#17934D", "#1A9850", "#229C52", 
                         "#2A9F54", "#31A355", "#39A757", "#41AB59", "#48AE5B", "#50B25D", 
                         "#58B65F", "#5FBA61", "#67BD63", "#6DC063", "#74C364", "#7AC665", 
                         "#81C865", "#87CB66", "#8ECE67", "#94D168", "#9BD468", "#A1D769", 
                         "#A7D96A", "#ACDB6E", "#B1DE71", "#B7E075", "#BCE278", "#C1E47B", 
                         "#C6E77E", "#CBE982", "#D0EB85", "#D5ED89", "#DAEF8D", "#DEF192", 
                         "#E2F297", "#E6F49C", "#E9F6A2", "#EDF7A7", "#F1F9AC", "#F5FAB1", 
                         "#F9FCB7", "#FDFEBC", "#FEFDBC", "#FEFAB7", "#FEF7B1", "#FEF4AC", 
                         "#FEF0A7", "#FEEDA2", "#FEEA9C", "#FEE797", "#FEE492", "#FEE18D", 
                         "#FDDC88", "#FDD784", "#FDD27F", "#FDCD7B", "#FDC877", "#FDC373", 
                         "#FDBE6F", "#FDB96A", "#FDB466", "#FDAF62", "#FCA95E", "#FBA25B", 
                         "#FA9C58", "#F99555", "#F88F52", "#F7884F", "#F6824C", "#F67B49", 
                         "#F57446", "#F46E43", "#F16840", "#EE613D", "#EB5B3B", "#E85538", 
                         "#E54F35", "#E34932", "#E0432F", "#DD3C2C", "#DA362A", "#D73027", 
                         "#D22B26", "#CD2626", "#C82126", "#C31D26", "#BE1826", "#B91326", 
                         "#B40E26", "#AF0926", "#AA0426", "#A50026")
  
  paleta_continuo_3 <- c("#440154FF", "#450558FF", "#46085CFF", "#470D60FF", "#471063FF", 
                         "#481467FF", "#481769FF", "#481B6DFF", "#481E70FF", "#482173FF", 
                         "#482576FF", "#482878FF", "#472C7AFF", "#472F7CFF", "#46327EFF", 
                         "#453581FF", "#453882FF", "#443B84FF", "#433E85FF", "#424186FF", 
                         "#404587FF", "#3F4788FF", "#3E4A89FF", "#3D4D8AFF", "#3C508BFF", 
                         "#3B528BFF", "#39558CFF", "#38598CFF", "#375B8DFF", "#355E8DFF", 
                         "#34608DFF", "#33638DFF", "#32658EFF", "#31688EFF", "#2F6B8EFF", 
                         "#2E6D8EFF", "#2D708EFF", "#2C718EFF", "#2B748EFF", "#2A768EFF", 
                         "#29798EFF", "#287C8EFF", "#277E8EFF", "#26818EFF", "#26828EFF", 
                         "#25858EFF", "#24878EFF", "#238A8DFF", "#228D8DFF", "#218F8DFF", 
                         "#20928CFF", "#20938CFF", "#1F968BFF", "#1F998AFF", "#1E9B8AFF", 
                         "#1F9E89FF", "#1FA088FF", "#1FA287FF", "#20A486FF", "#22A785FF", 
                         "#24AA83FF", "#25AC82FF", "#28AE80FF", "#2BB07FFF", "#2EB37CFF", 
                         "#31B67BFF", "#35B779FF", "#39BA76FF", "#3DBC74FF", "#41BE71FF", 
                         "#47C06FFF", "#4CC26CFF", "#51C56AFF", "#56C667FF", "#5BC863FF", 
                         "#61CA60FF", "#67CC5CFF", "#6DCD59FF", "#73D056FF", "#78D152FF", 
                         "#7FD34EFF", "#85D54AFF", "#8CD646FF", "#92D741FF", "#99D83DFF", 
                         "#A0DA39FF", "#A7DB35FF", "#ADDC30FF", "#B4DE2CFF", "#BBDE28FF", 
                         "#C2DF23FF", "#C9E020FF", "#D0E11CFF", "#D7E219FF", "#DDE318FF", 
                         "#E4E419FF", "#EBE51AFF", "#F1E51DFF", "#F7E620FF", "#FDE725FF")
  
  paleta_continuo_4 <- c("#000004FF", "#010107FF", "#02020BFF", "#030311FF", "#050417FF", 
                         "#07061CFF", "#090721FF", "#0C0926FF", "#0F0B2CFF", "#120D32FF", 
                         "#150E37FF", "#180F3EFF", "#1C1044FF", "#1F114AFF", "#221150FF", 
                         "#261257FF", "#2A115DFF", "#2F1163FF", "#331068FF", "#38106CFF", 
                         "#3C0F71FF", "#400F74FF", "#451077FF", "#491078FF", "#4E117BFF", 
                         "#51127CFF", "#56147DFF", "#5A167EFF", "#5D177FFF", "#611980FF", 
                         "#661A80FF", "#6A1C81FF", "#6D1D81FF", "#721F81FF", "#762181FF", 
                         "#792282FF", "#7D2482FF", "#822581FF", "#862781FF", "#8A2981FF", 
                         "#8E2A81FF", "#922B80FF", "#962C80FF", "#9B2E7FFF", "#9F2F7FFF", 
                         "#A3307EFF", "#A7317DFF", "#AB337CFF", "#AF357BFF", "#B3367AFF", 
                         "#B83779FF", "#BC3978FF", "#C03A76FF", "#C43C75FF", "#C83E73FF", 
                         "#CD4071FF", "#D0416FFF", "#D5446DFF", "#D8456CFF", "#DC4869FF", 
                         "#DF4B68FF", "#E34E65FF", "#E65163FF", "#E95562FF", "#EC5860FF", 
                         "#EE5C5EFF", "#F1605DFF", "#F2655CFF", "#F4695CFF", "#F66D5CFF", 
                         "#F7735CFF", "#F9785DFF", "#F97C5DFF", "#FA815FFF", "#FB8661FF", 
                         "#FC8A62FF", "#FC9065FF", "#FD9567FF", "#FD9A6AFF", "#FE9E6CFF", 
                         "#FEA36FFF", "#FEA873FF", "#FEAC76FF", "#FEB27AFF", "#FEB67DFF", 
                         "#FEBB81FF", "#FEC085FF", "#FEC488FF", "#FEC98DFF", "#FECD90FF", 
                         "#FED395FF", "#FED799FF", "#FDDC9EFF", "#FDE1A2FF", "#FDE5A7FF", 
                         "#FDEBABFF", "#FCEFB1FF", "#FCF4B6FF", "#FCF8BAFF", "#FCFDBFFF")
  
  paleta_continuo_5 <- c("#00FF00", "#00F900", "#00F400", "#00EF00", "#00EA00", "#00E500", 
                         "#00E000", "#00DA00", "#00D500", "#00D000", "#00CB00", "#00C600", 
                         "#00C100", "#00BC00", "#00B600", "#00B100", "#00AC00", "#00A700", 
                         "#00A200", "#009D00", "#009700", "#009200", "#008D00", "#008800", 
                         "#008300", "#007E00", "#007900", "#007300", "#006E00", "#006900", 
                         "#006400", "#005F00", "#005A00", "#005400", "#004F00", "#004A00", 
                         "#004500", "#004000", "#003B00", "#003600", "#003000", "#002B00", 
                         "#002600", "#002100", "#001C00", "#001700", "#001200", "#000C00", 
                         "#000700", "#000200", "#020000", "#070000", "#0C0000", "#120000", 
                         "#170000", "#1C0000", "#210000", "#260000", "#2B0000", "#300000", 
                         "#360000", "#3B0000", "#400000", "#450000", "#4A0000", "#4F0000", 
                         "#550000", "#5A0000", "#5F0000", "#640000", "#690000", "#6E0000", 
                         "#730000", "#790000", "#7E0000", "#830000", "#880000", "#8D0000", 
                         "#920000", "#970000", "#9D0000", "#A20000", "#A70000", "#AC0000", 
                         "#B10000", "#B60000", "#BC0000", "#C10000", "#C60000", "#CB0000", 
                         "#D00000", "#D50000", "#DA0000", "#E00000", "#E50000", "#EA0000", 
                         "#EF0000", "#F40000", "#F90000", "#FF0000")
  
  paleta_continuo_6 <- c("#FF0000FF", "#FF0F00FF", "#FF1F00FF", "#FF2E00FF", "#FF3D00FF", 
                         "#FF4D00FF", "#FF5C00FF", "#FF6B00FF", "#FF7A00FF", "#FF8A00FF", 
                         "#FF9900FF", "#FFA800FF", "#FFB800FF", "#FFC700FF", "#FFD600FF", 
                         "#FFE500FF", "#FFF500FF", "#FAFF00FF", "#EBFF00FF", "#DBFF00FF", 
                         "#CCFF00FF", "#BDFF00FF", "#ADFF00FF", "#9EFF00FF", "#8FFF00FF", 
                         "#80FF00FF", "#70FF00FF", "#61FF00FF", "#52FF00FF", "#42FF00FF", 
                         "#33FF00FF", "#24FF00FF", "#14FF00FF", "#05FF00FF", "#00FF0AFF", 
                         "#00FF1AFF", "#00FF29FF", "#00FF38FF", "#00FF47FF", "#00FF57FF", 
                         "#00FF66FF", "#00FF75FF", "#00FF85FF", "#00FF94FF", "#00FFA3FF", 
                         "#00FFB3FF", "#00FFC2FF", "#00FFD1FF", "#00FFE0FF", "#00FFF0FF", 
                         "#00FFFFFF", "#00F0FFFF", "#00E0FFFF", "#00D1FFFF", "#00C2FFFF", 
                         "#00B2FFFF", "#00A3FFFF", "#0094FFFF", "#0085FFFF", "#0075FFFF", 
                         "#0066FFFF", "#0057FFFF", "#0047FFFF", "#0038FFFF", "#0029FFFF", 
                         "#0019FFFF", "#000AFFFF", "#0500FFFF", "#1400FFFF", "#2400FFFF", 
                         "#3300FFFF", "#4200FFFF", "#5200FFFF", "#6100FFFF", "#7000FFFF", 
                         "#8000FFFF", "#8F00FFFF", "#9E00FFFF", "#AD00FFFF", "#BD00FFFF", 
                         "#CC00FFFF", "#DB00FFFF", "#EB00FFFF", "#FA00FFFF", "#FF00F5FF", 
                         "#FF00E6FF", "#FF00D6FF", "#FF00C7FF", "#FF00B8FF", "#FF00A8FF", 
                         "#FF0099FF", "#FF008AFF", "#FF007AFF", "#FF006BFF", "#FF005CFF", 
                         "#FF004CFF", "#FF003DFF", "#FF002EFF", "#FF001FFF", "#FF000FFF")
  
  #### Plot Continuous paletete ####
  listapaletas <- list(paleta_continuo_6, paleta_continuo_5, paleta_continuo_4, paleta_continuo_3, paleta_continuo_2,paleta_continuo_1,paleta_continuo)
  names(listapaletas) <- c("Colorful","RdBlGn","Magma","Viridis","RdYlGn", "Spectral", "Current")
  n <-sapply(listapaletas,function(x) length(x))
  nr <- length(listapaletas)
  nc <- max(n)
  
  ylim <- c(0,nr)
  oldpar <- graphics::par(mgp=c(2,0.25,0))
  on.exit(graphics::par(oldpar))
  plot(1,1,xlim=c(0,nc),ylim=ylim,type="n", axes=FALSE, bty="n",
       xlab="",ylab="",main = "pRocessomics color space - Continuous values")
  for(i in 1:nr){
    nj <- n[i]
    shadi <- listapaletas[[i]]
    graphics::rect(xleft=0:(nj-1), ybottom=i-1, xright=1:nj, ytop=i-0.2, col=shadi, border=shadi)
  }
  graphics::text(rep(-0.3,nr),(1:nr)-0.6, labels=names(listapaletas), xpd=TRUE, adj=1)
  
  palettes <- c("Spectral", "RdYlGn", "Viridis", "Magma", "RdBlGn","Colorful")
  choicetreat <- utils::select.list(choices = palettes,multiple = F,preselect = "Current",title = "\nPlease select a number of palette for heatmaps and other continuous scales[enter].")
  if(choicetreat=="Spectral") new_paleta_continuo <- paleta_continuo_1
  if(choicetreat=="RdYlGn") new_paleta_continuo <- paleta_continuo_2
  if(choicetreat=="Viridis") new_paleta_continuo <- paleta_continuo_3
  if(choicetreat=="Magma") new_paleta_continuo <- paleta_continuo_4
  if(choicetreat=="RdBlGn") new_paleta_continuo <- paleta_continuo_5
  if(choicetreat=="Colorful") new_paleta_continuo <- paleta_continuo_6
  revlvlcols <- utils::select.list(choices = truefalse,multiple = F, title = "\nDo you want to reverse the order of the colors?")
  if(revlvlcols=="Yes") new_paleta_continuo<-rev(new_paleta_continuo)
  
  new_continous_palette <- new_paleta_continuo
  class(new_continous_palette) <- c("color.pRo", "continous")
  new <- list(new_continous_palette)
  names(new) <- c("new_continuous_palette")
  list2env(new, envir = .GlobalEnv)
  
  #### Plot colors finales####
  graphics::par(mfrow=c(2,1),mar=c(0, 4.5, 1.5, 0))
  #upper plot
  listapaletasactual <- list(new_paleta_niveles, new_paleta_tratamientos)
  names(listapaletasactual) <- c("Omic Levels", "Treatments")
  n <-sapply(listapaletasactual,function(x) length(x))
  nr <- length(listapaletasactual)
  nc <- max(n)
  
  ylim <- c(0,nr)
  graphics::plot(1,1,xlim=c(0,nc),ylim=ylim,type="n", axes=FALSE, bty="n", xlab="",ylab="",main = "Current pRocessomics color space")
  for(i in 1:nr){
    nj <- n[i]
    shadi <- listapaletasactual[[i]]
    graphics::rect(xleft=0:(nj-1), ybottom=i-1, xright=1:nj, ytop=i-0.2, col=shadi,border="light grey")
  }
  graphics::text(rep(-0.2,nr),(1:nr)-0.6, labels=names(listapaletasactual), xpd=TRUE, adj=1)
  
  #lower plot
  listapaletasactual <- list(new_paleta_continuo)
  names(listapaletasactual) <- c("Continuous")
  n <-sapply(listapaletasactual,function(x) length(x))
  nr <- length(listapaletasactual)
  nc <- max(n)
  
  ylim <- c(0,nr)
  
  graphics::plot(1,1,xlim=c(0,nc),ylim=ylim,type="n", axes=FALSE, bty="n", xlab="",ylab="")
  for(i in 1:nr){
    nj <- n[i]
    shadi <- listapaletasactual[[i]]
    graphics::rect(xleft=0:(nj-1), ybottom=i-1, xright=1:nj, ytop=i-0.2, col=shadi,border=shadi)
  }
  graphics::text(rep(-0.2,nr),(1:nr)-0.6, labels=names(listapaletasactual), xpd=TRUE, adj=1)
  
  # #una vez que termina toda la funcion, la hacemos compatible con lo demas
  #  
  # new.aest <- paleta_tratamientos
  # class(new.aest) <- c("color.pRo", "treatment")
  # new.aesl <- paleta_niveles
  # class(new.aesl) <- c("color.pRo", "levels")
  # new.aesc <- paleta_continuo
  # class(new.aesc) <- c("color.pRo", "continuo")
  # nucolorspace <-list(new.aesc,new.aesl,new.aest)
  # names(nucolorspace) <- c("continuous_palette", "level_palette","treatment_palette")
  # list2env(nucolorspace)
  
  
}