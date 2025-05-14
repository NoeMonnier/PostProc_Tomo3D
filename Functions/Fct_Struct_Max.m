function  Im_Struct_Max = Fct_Struct_Max(Image)

% Fonction qui permet de conserver uniquement la structure d'aire max d'une image binaire

[sizey, sizex] = size(Image);
Im_Struct_Max = zeros(sizey, sizex);

ImageLabel = bwlabel(Image,4);
GrainData = regionprops(ImageLabel,'centroid','area');

for m = 1:length(GrainData)
    Structure(m) = GrainData(m).Area;
end
ind_max = find(Structure(:)==max(Structure(:))); % indice of largest structure
Im_Struct_Max(find(ImageLabel(:,:)==ind_max)) = 1;
