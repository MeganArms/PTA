classdef PTA < handle
   properties
       filename
       Option
       Molecule
       Frame
       Result
       Intensity
       Levy
       Surface
       Bulk
       Adsorption
       Desorption
       Aggregate
       Single
       Unknown
   end
   
   methods
       function obj = PTA(filename)
           if nargin >0
               obj.filename = filename;
           end
           obj.Option = struct;
       end
       
       function obj = set.Option(obj,opt)
           obj.Option = opt;
       end
       
       function setoption(obj)
           opt = struct;
           prompt = {'Spot radius (pixels)','Pixel size (nm)',...
               'Exclude region','Include only region',...
               'Connect distance threshold (nm)','Downsampling rate'...
               'Illumination correction (on or off)','Background'};
           dlg_title = 'Set the options';
           def = {'5','400','0','0','1200','1','off','1000'};
           answer = inputdlg(prompt,dlg_title,1,def);
           opt.spotR = str2double(answer{1});
           opt.pixelSize = str2double(answer{2});
           exclude = str2num(answer{3});
           if isempty(exclude)
               opt.exclude = false;
           else
               opt.exclude = exclude;
           end
           include = str2num(answer{4});
           if isempty(include)
               opt.include = false;
           else
               opt.include = include;
           end
           opt.connectDistance = str2double(answer{5});
           opt.ds = str2double(answer{6});
           opt.illumination = answer{7};
           opt.bg = str2double(answer{8});
           obj.Option = opt;
       end
   end
   
end
