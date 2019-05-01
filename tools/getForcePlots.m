clc; clear; clf;

% Get filenames
files=dir('force_????.txt');
nfiles=size(files,1);
filenames=repmat('force_????.txt',nfiles,1);  % Initialise

hold on;

for i=1:nfiles 
  filenames(i,:)=files(i).name;  % Extract filenames
end
legendnames=filenames(:,7:10);

for i=1:nfiles
  filename=filenames(i,:);
  data=dlmread(filename);  % Get data
  data(1,:)=[];  % Remove headers
  plot(data(:,1),data(:,2));
end

grid on;
legend(legendnames)
