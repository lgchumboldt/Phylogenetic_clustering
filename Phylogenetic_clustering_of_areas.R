###load packages###
library(maptools)
library(raster)
library(rgdal)
library(picante)
###Set working directory###
setwd("C:\\Users\\GIC 66\\Documents\\Andrea\\filodiversidad\\Beta_filodiversidad\\distribuciones_aves")
###Load phylogeny###
aves_iucn<-read.nexus("C:\\Users\\GIC 66\\Documents\\Andrea\\filodiversidad\\Aves\\Filogenias_aves_1590\\1k\\consenso_1k_aves_IUCN.nex")

#Folder containing distribution maps
distribution_maps_folder<-"C:\\Users\\GIC 66\\Documents\\Andrea\\filodiversidad\\Beta_filodiversidad\\distribuciones_aves"

distribution_files<-list.files(path=distribution_maps_folder, pattern= "*.shp$")
species_names<-sub(".shp","",distribution_files)
tabla<-as.data.frame(species_names)
colnames(tabla)<-"Grilla"

####Load mask file to be used in the run
colombia_shape<-"C:\\Users\\GIC 66\\Documents\\Andrea\\SIG\\Colombia\\adm\\COL_adm0.shp"
mascara<-readOGR(colombia_shape)

#Determine working resolution in degrees 
resolucion<-0.1


rasterize_species= function (x,mask=mascara) {
  r<-raster(ncol=1462,nrow=624)
  res(r)<-resolucion #resolution
  r<-crop(r, extent(mascara))
  values(r)<-0
  map<-readOGR(dsn=distribution_maps_folder,layer=x)
  r<-rasterize(map,r,1,update=T,background=0)
  r<-mask(r,mascara)
  valor<-unique(getValues(r))
  
  if(length(valor)==1&&is.na(valor)==TRUE){
    
    
  }
  else if (length(valor)==2&&valor[2]==0){
    
  }
  else {
    writeRaster(r,paste(x,".asc",sep=""))
    return (raster(paste(x,".asc",sep="")))
  }
}





layers<-lapply(species_names,rasterize_species)
names(layers)<-as.vector(tabla$Grilla)
layers[sapply(layers,is.null)]<-NULL
Stack_grilla<-stack(layers)



r<-raster(ncol=18,nrow=25)
extent(r)=extent(Stack_grilla)
r<-crop(r,extent(Stack_grilla))
res(r)<-res(Stack_grilla)
grilla=r
names(grilla)="grilla"
grilla[1:ncell(grilla)]<-1:ncell(grilla)

Stack_grilla<-stack(grilla,Stack_grilla)
marco<-as.data.frame(Stack_grilla)
marco<-na.omit(marco)
write.table(marco,"comunidades_aves_0_1.txt")

###########################################################################################
####Generation of evolutionary informed groups according to similarity between pixels######
###########################################################################################


#metodos=c("average","single","complete","ward","weighted")

marco1<-marco[,2:length(marco)]
rownames(marco1)<-marco[,1]
####Generate a distance matrix between all pixels. This can be done using several different functions: phylosor, unifrac seacrh for more
beta_PD_pixeles<-phylosor(marco1,aves_iucn)
####Perform a cluster analysis and plot a similitude dendrogram between pixels. Clustering can be donde using one of several methods: "average","single","complete","ward.D","weighted"
CLUSTER1=hclust(1-beta_PD_pixeles,  method ="ward")
plot(CLUSTER1)

##### Selecting the number of groups to use in the determination of evolutionary units (k).
groups<-20
###Plotting the selected groups on the dendrogram
rect.hclust(CLUSTER1,k=groups,border="blue")
###Plotting groups in a map
group_ras<-grilla
values(group_ras)<-NA
part_dendrogram<-cutree(CLUSTER1, k = gropus, h = NULL)
group_ras[as.numeric(names(part_dendrogram))]<-prueba3[names(part_dendrogram)]
plot(group_ras)
 


