function topofix(varargin)
%function topofix(varargin)
%
%------------------------------------------------------------------------
%GUI which can be used to hand-pick points and regions that are necessary
%to correct - original file is assumed to be from the confgrid routine.
%------------------------------------------------------------------------
%
% GUI is usually started in this way:
% topofix('size',[400 600])
%    The routine will now open the regional.depth/regional.grid files.
%    If regional.depth files cant be found, the routine will look for
%    the depthsaaaXbbb.uf file, where aaa and bbb are the grid sizes.
%    
%
% GUI can also be started this way
% topofix('depths',depth_matrix)
%    Where "depth_matrix" is a matrix specified in matlab - used as starting
%    point to the bathymetry.
%
%------------------------------------------------------------------------
% There are several options possible in the GUI:
% -Click on a point to select it (and probe it for values)
% -Click and drag to select a region
% -You can modify values of a point or region by setting a depth value in 
%  modify selection  and clicking "Modify"
% -You can change color axis of the plot
% -You can toggle the grid on and off
% -You can use the toolbar (ex zoom) - but switch it off before selecting new region
% -You can dump the bathymetry to file - this creates a topofix.depthsXXXxXXX.uf
%  and a topofix.regional.depth.[ab] files.
% -You can also dump the bathymetry to the matlab workspace.
%
% The last option makes it possible to do bathymetry processing in matlab
% without the guidance of topofix. Bathymetry smoothing, for instance
%------------------------------------------------------------------------
% Feel free to modify/improve in any way. User Interface 'hardening' may 
% be a place to start..
% 27th. Feb 2008, Knut Arild Liseter
%------------------------------------------------------------------------

% Makes the function "remember" these values
persistent origdepths lastdepths depths lon lat;

hasdepths=1==0;
infile=[];
idm=[];
jdm=[];

if nargin == 0
   %clear origdepths lastdepths depths lon lat;
   T=help('topofix');
   disp(T);

   disp('');
   disp('');
   disp('ERROR: topofix must be initialized with a grid size or a depth matrix!!')
   return;

% This is the case for callbacks
elseif nargin == 1 

   selector=varargin{1};

   % If
   if (isnumeric(selector))
	  if (selector==0)
		 disp('topofix must be initialized with a grid size or a depth matrix')
		 return
	  end
   else
	  disp('topofix must be initialized with a grid size or a depth matrix')
	  return
   end

else

   % This is an initialization statement
   selector=0;

	%  parse args pairvise
	iarg=1;
	nvarg=prod(size(varargin));
	while (iarg <= nvarg-1)
	   %varargin(iarg)
	   if (strcmp(varargin(iarg),'depths'))
		  hasdepths=1==1;
		  depths=varargin{iarg+1};
	   elseif (strcmp(varargin(iarg),'size'))
		  idm=varargin{iarg+1}(1);
		  jdm=varargin{iarg+1}(2);
	   else
		  disp(['uknown argument ' varargin{iarg}]);
		  return
	   end
	   iarg=iarg+2;
	end

	if (isempty(idm) & isempty(jdm) & isempty(depths))
	   disp('No appropriate arguments supplied')
	   return;
	elseif (~isempty(idm) & ~isempty(jdm) & isempty(depths))
	   disp(['Grid size supplied - will look for depth file with size ' num2str(idm) 'x' num2str(jdm)])
	elseif (isempty(idm) & isempty(jdm) & ~isempty(depths))
	   idm=size(depths,1);
	   jdm=size(depths,2);
	   hasdepths=1==1;
	   disp(['Depth matrix supplied - depth matrix is of size ' num2str(idm) 'x' num2str(jdm) ]);
	elseif (~isempty(idm) & ~isempty(jdm) & ~isempty(depths))
	   %Check that size matches
	   if (idm~=size(depths,1) | jdm~=size(depths,2))
		  disp('matrix and grid size mismatch - try calling with either depth matrix or grid sizes');
		  return
	   end
	else
	   disp(['Insufficient args supplied'])
	   return
	end
end % End argument parsing

if (nargin>1) 
T=help('topofix');
uiwait(msgbox(T,'modal'));
disp(T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch selector
case 0 % Initialize GUI application

	% Read grid and depths
	if (~hasdepths)
	   %depths=loada('regional.depth.a',1,idm,jdm);
      depths = readbathy(idm,jdm);
	   %lon=loada('regional.grid.a',1,idm,jdm);
	   %lat=loada('regional.grid.a',2,idm,jdm);
      obj=abfile('regional.grid.a','regional_grid');
	   lon=obj.getfield('plon',[],[]);
	   lat=obj.getfield('plat',[],[]);
	   if (isempty(depths))
        error('depths is empty - you are probably missing the regional.depth.[ab] file')
	   end
	   if (isempty(lon) | isempty(lat))
		  error('lon or lat is empty - you are probably missing the regional.grid.[ab] file')
	   end
      disp('Read regional.grid.[ab] and regional.depth.[ab]');
	end

	disp('Initializing GUI')
	fig = figure('WindowButtonUpFcn','topofix(4)','WindowButtonDownFcn','topofix(3)','ToolBar','figure') ; 


	% Text edit areas
	info.textxup  = uicontrol('Style','text','Position',[80 40 60 20]);
	uicontrol('Style','text','Position',[1 40 78 20],'String','x corner 1');
	info.textyup  = uicontrol('Style','text','Position',[80 20 60 20]);
	uicontrol('Style','text','Position',[1 20 78 20],'String','y corner 1');
	info.textxdwn = uicontrol('Style','text','Position',[80 80 60 20]);
	uicontrol('Style','text','Position',[1 80 78 20],'String','x corner 2');
	info.textydwn = uicontrol('Style','text','Position',[80 60 60 20]);
	uicontrol('Style','text','Position',[1 60 78 20],'String','y corner 2');


	% Depth, lon and lat probing
	info.depthprobe  = uicontrol('Style','text','Position',[80 140 60 20]);
	uicontrol('Style','text','Position',[1 140 78 20],'String','Depth:');
	info.lonprobe    = uicontrol('Style','text','Position',[80 120 60 20]);
	uicontrol('Style','text','Position',[1 120 78 20],'String','Longitude:');
	info.latprobe    = uicontrol('Style','text','Position',[80 100 60 20]);
	uicontrol('Style','text','Position',[1 100 78 20],'String','Latitude:');

	% Modify selection button
	uicontrol('Style','text','Position',[1 200 128 20],'String','Modify selection:');
	uicontrol('Style','pushbutton','Position',[80 180 50 20],'String','Modify!', ...
	          'Callback','topofix(5)');
	info.modvalue  = uicontrol('Style','edit','Position',[1 180 78 20]);
	% Isolated point check
	uicontrol('Style','pushbutton','Position',[0 320 120 20],'String','Check Isolated', ...
	          'Callback','topofix(13)');
	% Workspace Dump button
	uicontrol('Style','pushbutton','Position',[0 300 120 20],'String','Workspace Dump', ...
	          'Callback','topofix(11)');
	% File Dump button
	uicontrol('Style','pushbutton','Position',[0 280 120 20],'String','File Dump', ...
	          'Callback','topofix(9)');
	% Toggle grid button
	uicontrol('Style','pushbutton','Position',[0 380 120 20],'String','Toggle Grid', ...
	          'Callback','topofix(6)');
	% Reset view button
	uicontrol('Style','pushbutton','Position',[0 400 120 20],'String','Reset View', ...
	          'Callback','topofix(12)');
	% Restart button
	uicontrol('Style','pushbutton','Position',[0 260 120 20],'String','Reset depths', ...
	          'Callback','topofix(7)');
	% Undo
	uicontrol('Style','pushbutton','Position',[0 240 120 20],'String','Undo last', ...
	                           'Callback','topofix(8)');
    % Color axis
	uicontrol('Style','pushbutton','Position',[0 360 120 20],'String','Change Color axis', ...
	          'Callback','topofix(10)');
	info.clim1= uicontrol('Style','edit','Position',[0 340 60 20]);
	info.clim2= uicontrol('Style','edit','Position',[60 340 60 20]);

	% Init info on last button up/down - this info is also available from editbox
	info.xlastdwn=[];
	info.ylastdwn=[];
	info.xlastup=[];
	info.ylastup=[];

	info.selectpatch=[];
	info.selectpoint=[];



   % Set some variables and initialize map area
	origdepths=depths; % Keeps original depths so that you can start over
	lastdepths=depths; % Keeps last depths so that you can "undo"
   info.ax=axes('Position',[.3 .05 .65 .90]); hold on;
	I=find(depths==0); depths(I)=nan;
	info.P=pcolor(info.ax,depths');
	shading(info.ax,'flat'); 
	set(fig,'UserData',info, 'HandleVisibility','callback');
    
	
case 3 % WindowButtonwdownfcn
   fig=gcf;
 	info = get(fig,'UserData');
   [axx,axy]=  fig2axes(fig,depths) ;
	if ~isempty(axx) 
	   set(info.textxdwn,'String',num2str(axx));
	   set(info.textydwn,'String',num2str(axy));
	   info.xlastdwn=axx;
	   info.ylastdwn=axy;
	end 
	set(fig,'Userdata',info);

case 4 % WindowButtonwUpfcn
   %disp('hei down');
   fig=gcf;
 	info = get(fig,'UserData');
   [axx,axy]=  fig2axes(fig,depths) ;
	if ~isempty(axx) 
	   set(info.textxup,'String',num2str(axx));
	   set(info.textyup,'String',num2str(axy));
	   info.xlastup=axx;
	   info.ylastup=axy;
	end 

	% Probe if down and up points are the same - plot point if they are
	info=drawpoint(info,depths,lon,lat);

	%Mark selection with patch if down and up points differ and are non-empty
	info=drawpatch(info);
	set(fig,'Userdata',info);


case 5 % Modify button pressed
   fig=gcf;
   info = get(fig,'Userdata');

   % Check that current selection is ok
   if (~isempty(info.xlastdwn) & ~isempty(info.xlastup))
	  disp('Current Selection is OK!')
	  selectionx=min([info.xlastdwn info.xlastup]):max([info.xlastdwn info.xlastup]) ;
	  selectiony=min([info.ylastdwn info.ylastup]):max([info.ylastdwn info.ylastup]) ;
   else
     uiwait(msgbox('Nothing is selected'),'modal');
     disp('Nothing is selected');
	 return
   end 

   % Check that modification value is ok
   modval=str2num(get(info.modvalue,'String'));
   if (~isempty(modval))
     disp(['modification value is :' num2str(modval)])
   else
     disp('modification value is bad')
	 return
   end

   % Re-set depths to specifed values
   lastdepths=depths;
   depths(selectionx,selectiony)=modval;
   I=find(depths<.1); depths(I)=nan;

   % Re-init plot
   info=reinitaxes(info,depths,lon,lat);
   set(fig,'UserData',info);

case 6 % Toggle grid pressed
   fig=gcf;
   info = get(fig,'Userdata');
   edgecolor=get(info.P,'EdgeColor');
   if ( strcmp(edgecolor,'none'))
	  set(info.P,'Edgecolor','k');
	  set(info.P,'EdgeAlpha',.5);
   else
	  set(info.P,'Edgecolor','none');
   end

case 7 % Reset depths pressed
   % Clear topofix queue
   fig=gcf;
   info = get(fig,'Userdata');
   set(fig,'Userdata',info);
   buttonName=questdlg('This will reset the depths to what you started with. Is this what you want?', ...
                       'Reset Depths?','Yes','No','No');
   switch buttonName,
   case 'Yes'
      [depths,lastdepths]=startover(origdepths,depths,lon,lat);
   end



case 8 % Undo last  pressed
	fig=gcf;
	info = get(fig,'Userdata');

   % Recall last depth
   [depths,lastdepths]=startover(lastdepths,depths,lon,lat);

case 9 % file dump pressed
   fig=gcf;
   info = get(fig,'Userdata');
   I=find(isnan(depths)); depths(I)=0;
   writebathy(depths,'topofix');
   depths(I)=nan;

case 10 % caxis change pressed
   fig=gcf;
   info = get(fig,'Userdata');
   clim1=str2num(get(info.clim1,'String'));
   clim2=str2num(get(info.clim2,'String'));
   if (~isempty(clim1) & ~isempty(clim2))
	  cmax=max(max([clim1 clim2]));
	  cmin=min(min([clim1 clim2]));
	  set(info.ax,'CLim',[clim1 clim2])
   end 

case 11 % Workspace dump pressed
   assignin('base','depths_topofix',depths);
   disp('bathymetry dumped to workspace');
   uiwait(msgbox('bathymetry dumped to workspace as depths_topofix','modal'));

case 12 % Reset view pressed
   fig=gcf;
   info = get(fig,'Userdata');
   ax=info.ax;
   set(ax,'XLim',[0 size(depths,1)]);
   set(ax,'YLim',[0 size(depths,2)]);

case 13 % Check isolated points
   cnt=0;
   frst=1;
   found=1;
   iter=0;
   while  (found>0)
      found=0;
      newdepths=depths;
      for i=2:size(depths,1)-1
      for j=2:size(depths,2)-1
      if (~isnan(depths(i,j)))
         I=find(isnan(depths(i-1:i+1,j)));
         J=find(isnan(depths(i,j-1:j+1)));
         npts=prod(size(I))+prod(size(J));
         if (npts>=3) 
            disp(['Point ' num2str(i) ' ' num2str(j) ' has +3 land neighbours ']);
            newdepths(i,j)=nan;
            cnt=cnt+1;
            found=found+1;
         end 
      end
      end
      end
      depths=newdepths;
      iter=iter+1;
   end 

   % Needs axes2fig..
   if (cnt==0);
      disp('No isolated points found');
      uiwait(msgbox('No isolated points found','modal'));
   else
      fig=gcf;
      info = get(fig,'Userdata')
      disp(['Found and closed ' num2str(cnt) ' isolated points']);
      uiwait(msgbox(['Found and closed ' num2str(cnt) ' isolated points'],'modal'));;
      figure(fig);
      info=reinitaxes(info,depths,lon,lat);
      set(fig,'UserData',info);
   end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Auxillary functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [axx,axy]=  fig2axes(fig,depths)
	%Transform to axis coordinates
	currpt=get(fig,'Currentpoint');
	info = get(fig,'UserData');
	pos=get(fig,'Position');
	tmpunits=get(info.ax,'Units');
	set(info.ax,'Units','normalized');
	axpos=get(info.ax,'Position');
	XLim=get(info.ax,'XLim');
	YLim=get(info.ax,'YLim');
	set(info.ax,'Units',tmpunits);

	a=(XLim(2)-XLim(1))/(axpos(3)*pos(3)) ;
	b=XLim(1)-a*axpos(1)*pos(3);
	axx=a*currpt(1)+b;

	a=(YLim(2)-YLim(1))/(axpos(4)*pos(4)) ;
	b=YLim(1)-a*axpos(2)*pos(4);
	axy=a*currpt(2)+b;

	%size(depths);
	%axx=floor(axx)
	%axy=floor(axy)
	axx=max(1,floor(axx));
	axy=max(1,floor(axy));
	axx=min(axx,size(depths,1));
	axy=min(axy,size(depths,2));
   %size(depths)

	% Test should be changed
	try 
	   depths(floor(axx),floor(axy));
	catch
	   axx=[];
	   axy=[];
	   disp('Outside of domain')
	end

% Re-Init axes
function newinfo=reinitaxes(info,depths,lon,lat)
	% Keep axes limits
	XLim=get(info.ax,'XLim');
	YLim=get(info.ax,'YLim');
	cax=caxis(info.ax);
   edgecolor=get(info.P,'EdgeColor');
   edgealpha=get(info.P,'EdgeAlpha');


	delete(info.ax);
   info.ax=axes('Position',[.3 .05 .65 .90]); hold on;
	I=find(depths==0); depths(I)=nan;
	info.P=pcolor(info.ax,depths');
	set(info.ax,'XLim',XLim);
	set(info.ax,'YLim',YLim);
	shading(info.ax,'flat'); 
   caxis(info.ax,cax);
   set(info.P,'EdgeColor',edgecolor);
   set(info.P,'EdgeAlpha',edgealpha);

   % Re-draw selections
   info.selectpatch=[];
   info.selectpoint=[];
	info=drawpoint(info,depths,lon,lat);
	info=drawpatch(info);

	newinfo=info;

function newinfo=drawpatch(info)
	if (info.xlastdwn ~= info.xlastup & info.ylastdwn ~= info.ylastup & ~isempty(info.xlastup) ...
       & ~isempty( info.xlastdwn) )
	   axold=info.xlastdwn ;
	   ayold=info.ylastdwn ;
	   if (~isempty(info.selectpatch) )
		  delete(info.selectpatch)
		  info.selectpatch=[];
	   end
	   if (~isempty(info.selectpoint) )
		  delete(info.selectpoint)
		  info.selectpoint=[];
	   end

	   minx=min([axold info.xlastup]);
	   maxx=max([axold info.xlastup]);
	   miny=min([ayold info.ylastup]);
	   maxy=max([ayold info.ylastup]);

	   %P=patch([axold axold info.xlastup info.xlastup],[ayold info.ylastup info.ylastup ayold],[ .5 .5 .5]);
	   P=patch([minx minx maxx+1 maxx+1],[maxy+1 miny miny maxy+1],[ .5 .5 .5]);
	   set(P,'FaceAlpha',.5);
	   info.selectpatch=P;
	end
	newinfo=info;

function newinfo=drawpoint(info,depths,lon,lat)

	% Probe if down and up points are the same:
	if (info.xlastdwn == info.xlastup & info.ylastdwn == info.ylastup & ~isempty(info.xlastup) )
	   set(info.depthprobe,'String',num2str(depths(info.xlastup,info.ylastup),'%5.1f'));
	   set(info.lonprobe  ,'String',num2str(lon   (info.xlastup,info.ylastup),'%6.2f'));
	   set(info.latprobe  ,'String',num2str(lat   (info.xlastup,info.ylastup),'%6.2f'));

	   if (~isempty(info.selectpatch) )
		  delete(info.selectpatch);
		  info.selectpatch=[];
	   end

	   %mark selection with crosshair
	   if (~isempty(info.selectpoint) )
		  delete(info.selectpoint);
		  info.selectpoint=[];
	   end

	   P=plot(info.xlastup,info.ylastup,'m+');
	   set(P,'MarkerSize',20);
	   set(P,'LineWidth',2);
	   info.selectpoint=P;
	   %get(P)
	end
	newinfo=info;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [depths,lastdepths]=startover(newdepths,olddepths,lon,lat);
     fig=gcf;
	  info = get(fig,'UserData');
	  lastdepths=olddepths;
	  depths=newdepths;
	  info=reinitaxes(info,depths,lon,lat);
	  set(fig,'UserData',info);

