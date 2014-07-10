function A = xload(filename, display)

% Laden einer ASCII Datei
% Die ersten mit % gekennzeichneten Zeilen werden übersprungen bzw.
% mit A = xload(filename, 'v') angezeigt

comment = 1;
if nargin==1; display='off'; end %if    Keine Anzeige der Kommentare 
fid = fopen( filename, 'r');

if fid<0 error ([filename ' not found']); end

% Einlesen der Kommentarzeilen:
while comment;
 line = fgetl(fid);
 if length(line)  
   A    = sscanf(line, '%f')';
   rows = length(A);
   if line(1)=='%' & display=='v'; disp(line); 
   else if rows comment=0; end; % if rows
   end				% if line(1) == '%'
 end; 				% if length line
end;   				% while comment

% Einlesen der Daten:

A    = [ A; fscanf(fid, '%f',[rows,inf])']; 

fclose(fid);
