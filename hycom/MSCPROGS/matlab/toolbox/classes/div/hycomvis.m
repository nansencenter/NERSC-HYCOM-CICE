function hycomvis(infile,ftype)
% ----------------------------------------------------------------
% Matlab GUI to visualize hycom fields directly from hycom
% ([ab]) files. Use this GUI to read and visualize the fields 
% stored in hycom ".ab" files.
% It is called like this:
%
%    hycomvis(filename(s),filetype)
%
% Example:
%    hycomvis('archv.1991_228_00.a','archv')
%
% The input files can be any of four different types; "restart",
% "archv", "nersc_daily", or "nersc_weekly". If you supply 
% several % files they must be in a cellstring or a char array.
%
% ----------------------------------------------------------------
% The GUI is split into buttons, and a map. The map is simply 
% the model bathymetry. By clicking and dragging on the map you
% can select regions, sections or stations. 
%
% There are two modes for plotting - Section or Horizontal, 
% which can be set % by the radio buttons in the upper-left 
% corner. In Section mode you can only choose points or cross-
% sections. In Horizontal mode you can only choose regions. 
%
% The plot action buttons are enabled depending on the plot 
% type you select, and on the current selection in the map. 
% In addition there are menus for selecting variables to plot 
% and files to plot the variables from.
%
% You can produce several plots if you want. When you change the 
% active file, either through the next file button, or through 
% the file selection menu, the currently displayed plots will be
% updated to contain data from the newly selected file. To close 
% all plots from hycomvis, push the clear plots button.
%
% ----------------------------------------------------------------
% NB: The color axis variable is kept between plot updates. You can
%     use this to manually set the color axis in a section or field
%     plot figure, and plot updates will use this color axis. The 
%     same is true for the x and y limit properties for station plots
%
% ----------------------------------------------------------------
% hycomvis requires the "abfile" class to be installed. It quits 
% if it is not.
% Feel free to modify/improve in any way. 
% Created by Knut Arild Liseter, 17th Jan 2009. 
% ----------------------------------------------------------------


if (~exist('abfile'))
	disp('This program needs the abfile class - install it')
	return 
end

if (nargin==2)
  if (~strcmp(ftype,'nersc_daily') & ~strcmp(ftype,'nersc_weekly') ...
    & ~strcmp(ftype,'restart')  & ~strcmp(ftype,'archv') )
     disp('filetype must be one of nersc_daily , nersc_weekly, archv, or restart ');
     return;
  end
else
  disp('Usage - hycomvis(file, filetype)');
  return;
end

% Convert infiles to cell strings
if (ischar(infile) & min(min(size(infile)))>1 )
   'a'
	infile=cellstr(infile);
elseif (ischar(infile) & (size(infile,2)==1 | size(infile,1)==1 ))
   'b'
	% Hmm - this could be output from "ls" or something else 
	% - we parse it with strread and ' ' as delimiter - this returns cellstr
	infile=strread(infile,'%s');
	infile=sort(infile);
elseif (~iscellstr(infile))
	disp('input must be char array or cellstring array')
	return
end 


% These functions / tags are set by initgui 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Button/menu handles:
% initGUI           - GUI setup - placed at the very end of this file 
% wButtonUp         - Activated when mouse button released within axes
% wButtonDown       - Activated when mouse button pushed down within axes
% nextFile          - Next file Button, chooses next file in file list
% changeFile        - triggered when a new file is chosen
% buttonsCallback   - Handles button states when plot type is changed
% changeLev         - Handles level list when variable is changed
% plotFieldButton   - Handles actions when plot field button is clicked
% plotSectionButton - Handles actions when plot section button is clicked
% plotStationButton - Handles actions when plot station button is clicked
% clearButton       - Clears any open plots
% helpButton        - Opens instructions in a message box
%
% Tags - these are used to retrieve info on the GUI state:
% nextButton      - "Next File" button
% FilePopup       - file list
% sectionSwitch   - Tag for "section" plot type state
% clearButton     - Tag of the clear button
% helpButton      - Tag of the help button
% variablePopup   - Tag for variable selection menu 
% LevelPopup      - Tag for variable selection menu 
% fieldButton     - Tag for plot Field button
% stationButton   - Tag for plot Station button
% sectionButton   - Tag for plot Section button
% textXUp         - Tag for probe info - x position for last click release
% textYUp         - Tag for probe info - y position for last click release
% textXDown       - Tag for probe info - x position for last click push
% textYDown       - Tag for probe info - y position for last click push
% selectAxes      - Tag for depth / line / point / patch draw axes 
%
%Tags - these are used to retrieve info from the plot windows. ID
%is a unique identifier for a hycomvis session:
% HVFieldPlot+ID  - Tag denotes plot window is a field/horizontal plot
% HVStationPlot+ID- Tag denotes plot window is a station plot
% HVSectionPlot+ID- Tag denotes plot window is a section plot
%
%
%In addition - the GUI figure carries the following info in a structure 
%in  figure prop. "userdata"
%xlastdwn - last x pos of mouse push
%ylastdwn - last y pos of mouse push
%xlastup  - last x pos of mouse release
%ylastup  - last y pos of mouse release
%idm      - 1st data dim
%jdm      - 2nd data dim
%ftype    - file type
%String   - random string used to identify this hycomvis session
%


% Initialize GUI - at end of this file
initGUI(infile,ftype);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Button Selection handlers       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%% Handler detects when mouse button is released  %%%%%%%%%%%%%%%%%
% Sets main figure info for point, then draws points/patches/lines
function wButtonUp(source,eventdata,depths,lon,lat)
	[obj, fig] = gcbo; 
	updatePosition(fig,'buttonup',depths) ; % Update last position in userdata


	% Probe if down and up points are the same - plot point if they are. Only active 
	% if section plot type is selected
   if (get(findobj(fig,'Tag','sectionSwitch'),'Value')==1)
		if (drawpoint(fig,depths,lon,lat))
			% Enaple "plot station button" if drawpoint succeeded
			if (strcmp( get(findobj(fig,'Tag','stationButton'),'Enable'), 'off'))
			  set(findobj(fig,'Tag','stationButton'),'Enable','on');
			end
		end
	end

   % These actions are started by dragging. First is if section type plot 
	% is selected
   if (get(findobj(fig,'Tag','sectionSwitch'),'Value')==1)

		  %Mark selection with line if down and up points differ and are non-empty
		  if (drawline(fig))
			  % Enable section and station buttin if line draw succeeded
			  if (strcmp( get(findobj(fig,'Tag','sectionButton'),'Enable'), 'off'))
				 set(findobj(fig,'Tag','sectionButton'),'Enable','on');
			  end
			  if (strcmp( get(findobj(fig,'Tag','stationButton'),'Enable'), 'off'))
				 set(findobj(fig,'Tag','stationButton'),'Enable','on');
			  end
		  end

   % These actions are started by dragging. This section action is if field type plot 
	% is selected
   else 
		   %Mark selection with patch if down and up points differ and are non-empty
			if (drawpatch(fig))
				%Enable plot field button if this succeeded
				if (strcmp(get(findobj(fig,'Tag','fieldButton'),'Enable'),'off'))
					set(findobj(fig,'Tag','fieldButton'),'Enable','on');
				end
			else
				% Disable field button if we failed, also delete selectedPatch
				if (strcmp(get(findobj(fig,'Tag','fieldButton'),'Enable'),'on'))
					set(findobj(fig,'Tag','fieldButton'),'Enable','off');
				end
				clearselected(fig)
			end
   end




%%%%%%%%%%%%%%%%%% Handler detects when mouse button is pushed down %%%%%%%%%%%%%%%%%5
% Sets main figure info for point
function wButtonDown(source,eventdata,depths,lon,lat)
   [obj, fig] = gcbo; 
	updatePosition(fig,'buttondown',depths) ; % Update first position in userdata


%%%%%%%%%%%%%%%%%% Handler detects when plot field button is clicked %%%%%%%%%%%%%5
% Retrieves data and plots field using pcolor
function  plotFieldButton(source,eventdata,depths);
   [obj, fig] = gcbo; %Return calling object
 	info = get(fig,'UserData');
	filename = getActiveFile();
	fldname  = getActiveVar();
   level    = getActiveLevel();

   % Plot field
   pltfig=figure();
   plotField(pltfig,info,filename,fldname,level,depths) ;


%%%%%%%%%%%%%%%%%% Handler detects when plot section button is clicked %%%%%%%%%%%%%5
% Retrieves data and plots section using pcolor or line plot
function  plotSectionButton(source,eventdata,lon,lat);
   [obj, fig] = gcbo; %Return calling object
 	info = get(fig,'UserData');
	filename = getActiveFile();
	fldname = getActiveVar();

   % Plot section
   plotfig=figure();
   plotSection(plotfig,info,filename,fldname,lon,lat);


%%%%%%%%%%%%%%%%%% Handler detects when plot station button is clicked %%%%%%%%%%%%%5
% Retrieves data and plots station(s) using line plot
   function  plotStationButton(source,eventdata,lon,lat);
   [obj, fig] = gcbo; 
 	info = get(fig,'UserData');
	filename = getActiveFile();
   fldname = getActiveVar();

   % Plot section
   pltfig=figure();
   plotStations(pltfig,info,filename,fldname,lon,lat);


%%%%%% Handler detects when Mode is switched from field to section plot  %%%%%%%%%%%%
% Resets levels list (enabled or not) and clears anything that was selected
function buttonsCallback(source,eventdata)
   [obj, fig] = gcbo; 
   levlist=findobj(fig,'Tag','LevelPopup');
   if (get(findobj(fig,'Tag','sectionSwitch'),'Value')==1)
      set(levlist,'Enable','off'); % Switch level menu off for sections
   else
      set(levlist,'Enable','on');  % Switch level menu on for field plots
   end
   clearselected(fig); % Clear selection when context changes



%%%%%% Handler forces file change                         %%%%%%%%%%%%%%%%%%%%%%%%%%%
function nextFile(source,eventdata,lon,lat,depths)
	[obj, fig] = gcbo; 
   info=get(fig,'UserData');
   filelist=findobj(fig,'Tag','FilePopup');
	fldind=get(filelist,'Value');
	filenames=get(filelist,'String');
	fldind=mod(fldind,prod(size(filenames)))+1;
	set(filelist,'Value',fldind);

	% Trigger changeFile
	changeFile(source,eventdata);



%%%%%% Handler resets menus     when file     is changed %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Also updates any existing plots
function changeFile(source,eventdata)
	[obj, fig] = gcbo; 
   info=get(fig,'UserData');
	filename = getActiveFile()  ; 
	fldname  = getActiveVar()    ;
	level    = getActiveLevel()  ;
   disp(['File name changed ' filename])

	%Update variable list from contents of file
	obj = getFileObj(info.ftype,info.idm,info.jdm,fig);
	fldnames=getfieldnames(obj);
	varlist=findobj(fig,'Tag','VariablePopup'); 
	set(varlist,'String',fldnames);

   % Set variable to old list item if present
	I=find(strcmp(fldnames,fldname)==1);
	if (~isempty(I))
	   set(varlist,'Value',I);
	else
	   set(varlist,'Value',1);
	end

	% Trigger "change variable" function
	cv=changeVariable(source,eventdata)   ;

	filename = getActiveFile()  ; 

   % Update existing section plots
   splot=findobj('Tag',['HVSectionPlot.' info.String]);
   for i=1:prod(size(splot))
      plotSection(splot(i),info,filename);
   end

   % Update existing station plots
   splot=findobj('Tag',[ 'HVStationPlot.' info.String]);
   for i=1:prod(size(splot))
      plotStations(splot(i),info,filename);
   end

   % Update existing field plots
   splot=findobj('Tag',['HVFieldPlot.' info.String ]);
   for i=1:prod(size(splot))
      plotField(splot(i),info,filename) ;
   end


%%%%%% Handler resets level menu when variable is changed %%%%%%%%%%%%%%%%%%%%%%%%%%%
function success=changeVariable(source,eventdata)
   success=false;
   [obj, fig] = gcbo; info=get(fig,'UserData');
	filename = getActiveFile(); % Retrieved from GUI menu
	fldname  = getActiveVar();  % Retrieved from GUI menu

	% Get levels for  variable fldname in filemenu
   if (max(strcmp(info.ftype,{'nersc_daily' 'nersc_weekly' 'restart' 'archv'}))==1)
      obj=abfile(filename,info.ftype);
   end 
   flevels=getlevels(obj,fldname);


	% Set levels - if it is present in flevels return success...
	levlist=findobj(fig,'Tag','LevelPopup');
	oldlev=str2num(get(levlist,'String'));
	try 
		oldlev=reshape(oldlev,size(flevels));
		if (all(oldlev==flevels)) ; success=true ; end;
	end

	% Set Levels  and selected Value
   set(levlist,'String',num2str(flevels'));
	if (~success) ; set(levlist,'Value',1); end;



%%%%%%%%%% Clear any plots open in userdata of figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearButton(source,eventdata)
   [obj, fig] = gcbo; info=get(fig,'UserData');
   % Update existing section plots
   splot=findobj('Tag',['HVSectionPlot.' info.String]);
   for i=1:prod(size(splot))
      close(splot(i));
   end

   % Update existing station plots
   splot=findobj('Tag',[ 'HVStationPlot.' info.String]);
   for i=1:prod(size(splot))
      close(splot(i));
   end

   % Update existing field plots
   splot=findobj('Tag',['HVFieldPlot.' info.String ]);
   for i=1:prod(size(splot))
      close(splot(i));
   end

%%%%%%%%%% Print hycomvis help message %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function helpButton(source,eventdata)
   T=help('hycomvis');
   uiwait(msgbox(T,'modal'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Menu/button inquiry functions %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns active file in file menu
function filename = getActiveFile(fig);
   if (nargin==0)
		[obj, fig] = gcbo; 
	end
   filelist=findobj(fig,'Tag','FilePopup');
	fldind=get(filelist,'Value');
	filenames=get(filelist,'String');
   filename=filenames{fldind};

% Returns active var in var menu
function varname = getActiveVar(fig);
   if (nargin==0)
		[obj, fig] = gcbo; 
	end
   varlist=findobj(fig,'Tag','VariablePopup');
	varind=get(varlist,'Value');
	varnames=get(varlist,'String');
   varname=char(varnames{varind});


% Returns active level in level menu
function level = getActiveLevel();
   [obj, fig] = gcbo; 
   levlist=findobj(fig,'Tag','LevelPopup');
	levind=get(levlist,'Value');
	levels=get(levlist,'String');
   level=str2num(levels(levind,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Functions handling selection / drawing%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draws a patch which is kept in the info field of figure 
% Also modifies enabled status of buttons
function success=drawpatch(fig)
   info=get(fig,'UserData');
	success=false;
	if ( info.xlastdwn ~= info.xlastup & info.ylastdwn ~= info.ylastup & ...
        ~isempty(info.xlastup)  & ~isempty( info.xlastdwn) )
	   clearselected(fig);
	   minx=min([info.xlastdwn info.xlastup]);
	   maxx=max([info.xlastdwn info.xlastup]);
	   miny=min([info.ylastdwn info.ylastup]);
	   maxy=max([info.ylastdwn info.ylastup]);
	   P=patch([minx minx maxx+1 maxx+1],[maxy+1 miny miny maxy+1],[ .5 .5 .5], ...
              'Tag','selectedPatch');
	   set(P,'FaceAlpha',.5);
		success=true;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draws a line which is kept in the info field of mainfigure
function success=drawline(fig)
   info=get(fig,'UserData');
	success=false;
   if (info.xlastdwn ~= info.xlastup & info.ylastdwn ~= info.ylastup &  ...
       ~isempty(info.xlastup)  & ~isempty( info.xlastdwn) )
	  minx=info.xlastdwn;
	  maxx=info.xlastup;
	  miny=info.ylastdwn;
	  maxy=info.ylastup;
	  clearselected(fig);
	  L=line([info.xlastdwn info.xlastup],[info.ylastdwn info.ylastup], ...
            'Color','m','LineWidth',2,'Tag','selectedLine');
	  success=true;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draws a point which is kept in the info field of main figure (switch to tag?)
% Also updates the "probing" values
function success=drawpoint(fig,depths,lon,lat)
   success=false;
   info=get(fig,'UserData');
	% Probe if down and up points are the same:
	if (info.xlastdwn == info.xlastup & info.ylastdwn == info.ylastup & ...
       ~isempty(info.xlastup) )
	   set(findobj(fig,'Tag','depthProbe'),'String',...
          num2str(depths(info.xlastup,info.ylastup),'%5.1f'));
	   set(findobj(fig,'Tag','lonProbe'),'String',...
          num2str(lon(info.xlastup,info.ylastup)));
	   set(findobj(fig,'Tag','latProbe'),'String',...
          num2str(lat(info.xlastup,info.ylastup)));
	   clearselected(fig);
	   P=plot(info.xlastup,info.ylastup,'m+','Tag','selectedPoint');
	   set(P,'MarkerSize',20);
	   set(P,'LineWidth',2);
		success=true;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function clears all selected objects kept in the userdata part of the figure and
% removes them from the map. Also disables buttons, since we no longer have a selection
 function clearselected(fig);
   selectpatch=findobj(fig,'Tag','selectedPatch');
   if(~isempty(selectpatch)) 
      delete(selectpatch);
   end
	if (~strcmp( get(findobj(fig,'Tag','fieldButton'),'Enable'), 'off'))
      set(findobj(fig,'Tag','fieldButton'),'Enable','off');
	end
   selectpoint=findobj(fig,'Tag','selectedPoint');
	if (~isempty(selectpoint) )
	   delete(selectpoint)
	end
   if (~strcmp( get(findobj(fig,'Tag','stationButton'),'Enable'), 'off'))
     set(findobj(fig,'Tag','stationButton'),'Enable','off');
   end
   selectline=findobj(fig,'Tag','selectedLine');
	if (~isempty(selectline) )
	   delete(selectline)
	end
   if (~strcmp( get(findobj(fig,'Tag','sectionButton'),'Enable'), 'off'))
     set(findobj(fig,'Tag','sectionButton'),'Enable','off');
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Auxillary functions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Maps figure point to depth matrix point  - uses currentpoint fig prop %%%%%%%%
function [axx,axy]=  fig2axes(fig)
	%Transform to axis coordinates
	currpt=get(fig,'Currentpoint');
	info = get(fig,'UserData'); 
	pos=get(fig,'Position');

   ax=findobj(fig,'Tag',[ 'selectAxes.' info.String]);
	tmpunits=get(ax,'Units');
	set(ax,'Units','normalized');
	axpos=get(ax,'Position');
	XLim=get(ax,'XLim');
	YLim=get(ax,'YLim');
	set(ax,'Units',tmpunits);

	a=(XLim(2)-XLim(1))/(axpos(3)*pos(3)) ;
	b=XLim(1)-a*axpos(1)*pos(3);
	axx=a*currpt(1)+b;

	a=(YLim(2)-YLim(1))/(axpos(4)*pos(4)) ;
	b=YLim(1)-a*axpos(2)*pos(4);
	axy=a*currpt(2)+b;

	axx=max(1,min(floor(axx),info.idm));
	axy=max(1,min(floor(axy),info.jdm));



%%%%%%%%%% Update map position in userdata of figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatePosition(fig,direction,depths)
	info = get(fig,'UserData');
   [axx,axy]=  fig2axes(fig) ; % Get point on axes
   if ~isempty(axx) 
		if (strcmp(direction,'buttondown'))
			set(findobj(fig,'Tag','textXDown'),'String',num2str(axx));
			set(findobj(fig,'Tag','textYDown'),'String',num2str(axy));
			info.xlastdwn=axx;
			info.ylastdwn=axy;
		elseif (strcmp(direction,'buttonup'))
			set(findobj(fig,'Tag','textXUp'),'String',num2str(axx));
			set(findobj(fig,'Tag','textYUp'),'String',num2str(axy));
			info.xlastup=axx;
			info.ylastup=axy;
		else
			disp(['Unknown direction '  direction])
		end 
   end 
   set(fig,'Userdata',info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function returns data  points along a section - straight line on model grid
function [x,y,ind]=getsectionpoints(x1,y1,x2,y2,idm,jdm);
   maxdiff=max(abs(x1-x2),abs(y1-y2))+1;
	x=round(linspace(x1,x2,maxdiff));
	y=round(linspace(y1,y2,maxdiff));
   ind = sub2ind([idm jdm],x,y) ;
	[x,y]=ind2sub([idm jdm],ind);
	if (nargout==1)
	   x=ind;
	end
	   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Data retrieval routines             %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Function returns data values using data matrix points along a section (x,y)
function [fld,prs] = getsectiondata(x,y,infile,ftype,fldname,idm,jdm);
   % Retrieve fields from infile
   if (strcmp(ftype,'nersc_daily'))
      obj=abfile(infile,ftype);
      prs=getpoint(obj,x,y,'pres',[],[]);
      fld=getpoint(obj,x,y,fldname,[],[]); 
   elseif (strcmp(ftype,'nersc_weekly'))
      obj=abfile(infile,ftype);
      prs=getpoint(obj,x,y,'pres',[],[]);
      for k=2:size(prs,1)
         prs(k,:)=prs(k,:)+prs(k-1,:);
      end
      fld=getpoint(obj,x,y,fldname,[],[]);
   elseif (strcmp(ftype,'restart'))
      obj=abfile(infile,ftype);
      prs=getpoint(obj,x,y,'dp',[],1);
      for k=2:size(prs,1)
         prs(k,:)=prs(k,:)+prs(k-1,:);
      end
      fld=getpoint(obj,x,y,fldname,[],1); 
   elseif (strcmp(ftype,'archv'))
      obj=abfile(infile,ftype);
      prs=getpoint(obj,x,y,'thknss',[],1);
      for k=2:size(prs,1)
         prs(k,:)=prs(k,:)+prs(k-1,:);
      end
      fld=getpoint(obj,x,y,fldname,[],1); 
   end

   %TODO: Should fix this at "reader" end
   I=find(fld > 1e20); fld(I)=nan;
   I=find(prs > 1e20); prs(I)=0.;

   % Add a layer on top
   if (size(prs)==size(fld))
      prs=[zeros(1,size(prs,2));  prs];
      fld=[fld(1,:);  fld];
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function returns fields for chosen field and level
function fld = getfld(infile,ftype,fldname,vlevel,idm,jdm);
   % Retrieve fields from infile
   if (max(strcmp(ftype,{'nersc_daily' 'nersc_weekly' 'restart' 'archv'}))==1)
      obj=abfile(infile,ftype);
      fld=getfield(obj,fldname,vlevel,1);
   end



% Returns file object corr to current file selected
function obj = getFileObj(ftype,idm,jdm,fig);
   if (nargin==3)
		[obj, fig] = gcbo; 
	end
	filename = getActiveFile(fig);
	if (max(strcmp(ftype,{'restart' 'nersc_daily' 'nersc_weekly' 'archv'}))==1)
		obj=abfile(filename,ftype);
	else 
		disp(['Unknown file type ' ftype ]);
	end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Plotting routines                   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  plotStations(pltfig,info,filename,fldname,lon,lat);
% plots Stations 

   figure(pltfig);
   xlim=[];
   ylim=[];
   % First plot - everything calculated from scratch
   if (nargin==6)
      [x,y,ind]=getsectionpoints(info.xlastdwn,info.ylastdwn, ...
                           info.xlastup ,info.ylastup , ...
                           info.idm,info.jdm);
   % Update plot - use existing data
   elseif (nargin==3)
      myinfo=get(pltfig,'UserData');
      x      =myinfo.x;
      y      =myinfo.y;
      ind    =myinfo.ind;
      lon    =myinfo.lon;
      lat    =myinfo.lat;
      fldname=myinfo.fldname;
      obj=findobj(pltfig,'Tag',['HVStationPlot.', info.String]);
      if (~isempty(obj));
         xlim=get(obj,'XLim');
         ylim=get(obj,'YLim');
      end
   end
   [fld,prs] = getsectiondata(x,y,filename,info.ftype,fldname,info.idm,info.jdm);

   x2=repmat(1:prod(size(x)),size(prs,1),1);
   if (max(max(prs))>10.*9806)
	  dfac=9806;
   else
	  dfac=1;
   end
   if (size(fld)==size(prs)) ; 
		I=find(fld>1e20); fld(I)=nan;
		I=find(prs>1e20); prs(I)=nan;
      for k=1:size(prs,2)
         %prs(k,:)
         I=find(max(prs(:,k))-prs(:,k)<dfac);
         fld(I,k)=nan;
      end
      plot(fld,-prs,'.-'); shading flat;
      %set(gca,'FontSize',14);
      %set(gca,'FontWeight','bold');
      ylabel('Depth[m]'); 
   else
		I=find(fld>1e20); fld(I)=nan;
      plot(x2,fld,'LineWidth',2); 
      if (~isempty(xlim) & ~isempty(ylim))
         set(gca,'XLim',xlim);
         set(gca,'YLim',ylim);
      end
      %set(gca,'FontSize',14);
      %set(gca,'FontWeight','bold');
      XT=get(gca,'XTick');
      for i=1:prod(size(XT))
         if XT(i) < 1 | XT(i) > prod(size(ind));
            XTL{i}='';
         else
            XTL{i}=num2str(lon(ind(XT(i))),'%7.2f'  );
         end
      end
      set(gca,'XTickLabel',XTL); 
      xlabel('Longitude[degrees east]'); 
      grid on;
   end

   % Set figure info
   myinfo.x=x;
   myinfo.y=y;
   myinfo.ind=ind;
   myinfo.lon=lon;
   myinfo.lat=lat;
   myinfo.fldname=fldname;
   set(pltfig,'Tag',[ 'HVStationPlot.' info.String]);
   set(pltfig,'UserData',myinfo);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSection(pltfig,info,filename,fldname,lon,lat);
% Plots section - tags figure window containing plot and puts the figure handler
%into  info struct. 

   %Retrieve section data 
   figure(pltfig);

   % First plot - everything calculated from scratch
   if (nargin==6)
      [x,y,ind]=getsectionpoints(info.xlastdwn,info.ylastdwn, ...
                           info.xlastup ,info.ylastup , ...
                           info.idm,info.jdm);

      cax=[];
      xlim=[];
      ylim=[];
   % Update plot - use existing data
   elseif (nargin==3)
      myinfo=get(pltfig,'UserData');
      x      =myinfo.x;
      y      =myinfo.y;
      ind    =myinfo.ind;
      lon    =myinfo.lon;
      lat    =myinfo.lat;
      fldname=myinfo.fldname;
      cax=caxis;
      xlim=get(gca,'XLim');
      ylim=get(gca,'YLim');
   end

   % Retrieve field and depth levels (if applicable)
   [fld,prs] = getsectiondata(x,y,filename,info.ftype,fldname,info.idm,info.jdm);
	x2=repmat(1:prod(size(x)),size(prs,1),1);

   % Actual plot
   if (size(fld)==size(prs)) ; 
      dfac=1;
      if (max(max(prs))>10.*9806)
         dfac=1/9806;
      end
      P=pcolor(flipud(x2),-flipud(prs)*dfac,flipud(fld)); shading flat;
      obj=gca;
      if (~isempty(cax)); caxis(obj,cax); end
      if (~isempty(xlim)); set(obj,'XLim',xlim); end
      if (~isempty(ylim)); set(obj,'YLim',ylim); end
      colorbar;
		disp(['Caxis is ' num2str(caxis) ])
		disp(['XLim  is ' num2str(get(obj,'XLim')) ])
		disp(['YLim  is ' num2str(get(obj,'YLim')) ])
   else
		I=find(fld>1e26); fld(I)=nan;
      plot(x2,fld,'LineWidth',2);
      grid on;
   end

   %Labels (longitude for now)
   %set(gca,'FontSize',14);
   %set(gca,'FontWeight','bold');
   XT=get(gca,'XTick');
   for i=1:prod(size(XT))
      if XT(i) <= 0 | XT(i) > prod(size(ind));
         XTL{i}='';
      else
         XTL{i}=num2str(lon(ind(XT(i))),'%7.2f'  );
      end
   end
   set(gca,'XTickLabel',XTL); 
   xlabel('Longitude[degrees east]'); 
   if (size(fld)==size(prs)) ; 
      ylabel('Depth[m]'); 
   end

   % Set figure info
   myinfo.x=x;
   myinfo.y=y;
   myinfo.ind=ind;
   myinfo.lon=lon;
   myinfo.lat=lat;
   myinfo.fldname=fldname;
   set(pltfig,'Tag',['HVSectionPlot.' info.String]);
   set(pltfig,'UserData',myinfo);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotField(pltfig,info,filename,fldname,level,depths) ;
% Plots field - tags figure window containing plot and puts the figure handler
%into  info struct. 

	figure(pltfig); 
   if (nargin==6)
      minx=min([info.xlastdwn info.xlastup]);
      maxx=max([info.xlastdwn info.xlastup]);
      miny=min([info.ylastdwn info.ylastup]);
      maxy=max([info.ylastdwn info.ylastup]);
      cax=[];
   elseif (nargin==3)
      myinfo=get(pltfig,'UserData');
      minx   =myinfo.minx;
      maxx   =myinfo.maxx;
      miny   =myinfo.miny;
      maxy   =myinfo.maxy;
      level  =myinfo.level;
      fldname=myinfo.fldname;
      depths =myinfo.depths ;
      cax=caxis;
   end
   fld = getfld(filename,info.ftype,fldname,level,info.idm,info.jdm);
   I=find(depths<.1 | depths > 1e20 | isnan(depths)); fld(I)=nan;

   P=pcolor(minx:maxx,miny:maxy,fld(minx:maxx,miny:maxy)'); shading flat;
   if (~isempty(cax))
      caxis(cax);
   end
   colorbar;

   % Set figure info
   myinfo.minx  =minx;
   myinfo.maxx  =maxx;
   myinfo.miny  =miny;
   myinfo.maxy  =maxy;
   myinfo.depths=depths; % Not ideal to store large matrix in userdata
   myinfo.fldname=fldname;
   myinfo.level  =level  ;
   set(pltfig,'Tag',['HVFieldPlot.' info.String] );
   set(pltfig,'UserData',myinfo);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GUI initialization routine - sets tags and connects handles %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initGUI(infile,ftype);
	disp('Initializing GUI')

	% Get lon and lat
	rg=abfile('regional.grid.a','regional_grid');
	lon=getfield(rg,'plon',[],[]);
	lat=getfield(rg,'plat',[],[]);
	idm=size(lon,1);
	jdm=size(lon,2);
	rd=abfile('regional.depth.a','raw');
	depths=getfield(rd,[],1,[]);
	if (isempty(depths))
		disp('depths is empty - you are probably missing the regional.depth.a file')
	end
	if (isempty(lon) | isempty(lat))
		disp('lon or lat is empty - you are probably missing the regional.grid.a file')
	end

	fig = figure('WindowButtonUpFcn'  ,{@wButtonUp  ,depths,lon,lat}, ...
					 'WindowButtonDownFcn',{@wButtonDown,depths,lon,lat}, ...
					 'ToolBar','figure') ; 


	% Create a list box with files
	p1=uipanel('Title','Files','Position',[.01 .71 .30 .15 ]);
	uicontrol('Parent',p1,'Style','pushbutton', 'Units','normalized', ...
             'Position', [ .05 .05 .9 .5], 'String','Next File', ...
             'Callback', {@nextFile,lon,lat,depths}, ...
				 'Enable','on','Tag','nextButton');
	uicontrol('Parent',p1,'Style','popupmenu', 'Units','normalized',  ...
             'Position',[.05 .55 .9 .5 ], 'String',infile,  'Enable','on', ...
             'Tag','FilePopup',  'Callback', @changeFile );


	% Create the section/field  button group.
	p2=uipanel('Title','Plot Type','Position',[.01 .90 .50 .10 ]);
	h = uibuttongroup('visible','off','Units','normalized','Position',[0 0 1 1], ...
                      'Tag','buttonGroup', 'Parent',p2, 'BorderType','none');
	b0 = uicontrol('Style','Radio','String','Section/Station','Units','Normalized',...
						'pos',[.01 .01 .5 .9],'parent',h,'HandleVisibility','on', ...
                  'Tag','sectionSwitch');
	b1 = uicontrol('Style','Radio','String','Horizontal','Units','Normalized',...
						'pos',[.51 .01  .5 .9],'parent',h,'HandleVisibility','on');
	set(h,'SelectedObject',b0);  % No selection
	set(h,'SelectionChangeFcn',@buttonsCallback);
	set(h,'Visible','on');

   % Create the ``clear plots'' button
	uicontrol('Style','pushbutton', 'Units','Normalized', ...
             'Position',[0.60 0.90 .15 .08], 'String','Clear plots!', ...
             'Enable','on', 'Tag', 'clearButton', 'Callback',@clearButton);

   % Create the ``help'' button
	uicontrol('Style','pushbutton', 'Units','Normalized', ...
             'Position',[0.80 0.90 .15 .08], 'String','Help!',  ...
             'Enable','on', 'Tag', 'helpButton', 'Callback',@helpButton);

	% Create a list box with variables
	p3=uipanel('Title','Variable and Level Selection','Position',[.01 .51 .30 .15 ]);
	obj = getFileObj(ftype,idm,jdm,fig);
	fldnames=getfieldnames(obj);
	varlist    =uicontrol('Style','popupmenu','Units','Normalized', ...
                         'Position',[0.05 .55 .9 .45], 'String',fldnames, ...
                         'Enable','on','Callback',{@changeVariable},  ...
								 'Tag','VariablePopup','Parent',p3);  
	%Get variable name of initial entry of the uicontrol varlist "VariablePopup"
	fldname = getActiveVar(fig);
	flevels=getlevels(obj,fldname);
	uicontrol('Style','popupmenu','Units','Normalized', ...
             'Position',[0.05 .05 .9 .45],'String',flevels, ...
				 'Enable','off','Tag','LevelPopup','Parent',p3);


	% Section/Station/Field plot buttons
	p4=uipanel('Title','Plot Actions','Position',[.01 .27 .30 .20 ]);
	uicontrol('Style','pushbutton', 'Units','Normalized',...
             'Position',[0.05 0.66 .90 .33], 'String','Plot Field  ', ...
             'Callback',{@plotFieldButton,depths}, 'Enable','off', ...
				 'Tag', 'fieldButton','Parent',p4);
	uicontrol('Style','pushbutton','Units','Normalized',...
             'Position',[ 0.05 0.33 .90 .33 ], 'String','Plot Section',...
             'Callback',{@plotSectionButton,lon,lat},'Enable','off', ...
				 'Tag','sectionButton','Parent',p4);
	uicontrol('Style','pushbutton','Units','Normalized',...
             'Position',[0.05 0.00 .90 .33], 'String','Plot Station(s)',...
             'Callback', {@plotStationButton,lon,lat}, 'Enable','off', ...
				 'Tag','stationButton', 'Parent',p4);

	% Text edit areas
	p5=uipanel('Title','Location','Position',[.01 .01 .30 .20 ]);
	uicontrol('Style','text','Units','Normalized','Position',[.01 .2 .4 .2], ...
	          'HorizontalAlignment','right','String','Point 1:','Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.01 .0 .4 .2], ...
	          'HorizontalAlignment','right','String','Point 2:','Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.01 .8 .4 .2], ...
	          'HorizontalAlignment','right','String','Depth:','Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.01 .6 .4 .2], ...
	          'HorizontalAlignment','right','String','Longitude:','Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.01 .4 .4 .2], ...
	          'HorizontalAlignment','right','String','Latitude:','Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.41 .0 .3 .2 ], ...
	          'Tag','textXUp'   ,'Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.71 .0 .3 .2 ], ...
	          'Tag','textYUp'   ,'Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.41 .2 .3 .2 ], ...
	          'Tag','textXDown' ,'Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.71 .2 .3 .2 ], ...
	          'Tag','textYDown' ,'Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.41 .8 .6 .2 ], ...
	          'Tag','depthProbe','Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.41 .6 .6 .2 ], ...
	          'Tag','lonProbe'  ,'Parent',p5);
	uicontrol('Style','text','Units','Normalized','Position',[.41 .4 .6 .2 ], ...
	          'Tag','latProbe'  ,'Parent',p5);

   % Random string used to identify this main window (in case someone opens hycomvis for more 
   % than one file set / model)
   LetterStore = char(97:122); % string containing all allowable letters (in this case lower case only)
   Idx = randperm(length(LetterStore));
   String = LetterStore(Idx(1:6));


	% Init info on last button up/down 
	info.xlastdwn=[];
	info.ylastdwn=[];
	info.xlastup=[];
	info.ylastup=[];
	info.idm    =idm;
	info.jdm    =jdm;
	info.ftype=ftype;
	info.String=String;

	% Plot depths map
	ax=axes('Position',[.35 .05 .6  .80],'Tag',['selectAxes.' info.String] ); hold on;
	I=find(depths==0); depths(I)=nan;
	pcolor(ax,depths');
	shading(ax,'flat'); 
	set(fig,'UserData',info);
	set(fig,'HandleVisibility','callback');
	 
