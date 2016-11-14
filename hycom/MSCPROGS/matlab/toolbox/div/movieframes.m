function movieframes(varname,cax,lcolor,maskcond)
% Read selected variable from a tmp1.nc file (use m2nc), and create png-files 
% from it
warning off MATLAB:FINITE:obsoleteFunction;
global MAP_PROJECTION;
if (isempty(MAP_PROJECTION))
	disp('movieframes needs a call to m_proj first');
	return
	%Examples
	%m_proj('stereographic','lon',-5,'lat',79.5,'rad',7,'rec','on','rot',-13); - Close-up Fram Strait
	%m_proj('lambert','lon',[-20 20],'lat',[73.8 83],'rec','on'); - Close-up Fram Strait
end 


ncinfo=nc_info('tmp1.nc');
rdimlen=ncinfo.Dimension(3).Length;
nx=ncinfo.Dimension(1).Length;
ny=ncinfo.Dimension(2).Length;

lon  =nc_varget('tmp1.nc','longitude');
lat  =nc_varget('tmp1.nc','latitude');
depth=nc_varget('tmp1.nc','depth');
[x,y]=m_ll2xy(lon,lat);
%for i=1:rdimlen
frst=1==1;
for i=1:rdimlen
%for i=150:150
	fld=nc_varget('tmp1.nc',varname,[i-1 0 0],[1 -1 -1]);
	I=find(isnan(depth)); fld(I)=nan;

	if nargin==4
		I=eval([ 'find(fld' maskcond ');']);
		fld(I)=nan;
	end



	alphamask=zeros(size(fld));
	I=find(isnan(fld)) ; alphamask(I)=0;
	P=pcolor(x,y,fld); set(P,'EdgeColor','none')
	%set(P,'AlphaData',alphamask)
	%set(gca,'Color','none')
	%set(gcf,'Color','none')
	%set(gcf, 'InvertHardCopy', 'off');
	%get(P)
	caxis(cax);
	if (frst) 
		m_gshhs_i('patch',lcolor);
		m_gshhs_i('save','mycoast.mat');
		frst=0==1;
	else
		m_usercoast('mycoast.mat','patch',lcolor);
	end
	m_grid;
	cii=num2str(i,'%5.5d');
	fname=[ 'Pics/' varname '-' cii '.png' ];
	disp(fname)
	print('-dpng',fname);
end


