function [f] = splot(seq, mint, maxt)
%
% SEQUENCE.TOOLS.SPLOT
%
%	Plot the time trace for qc purpose of a generated sequence
%
% INPUT
%      seq  sequence structure generated with +sequence framework
%      mint [optional] minimum time range to show x-axis from in s
%      maxt [optional] maximum time range to show x-axis to in s
%
% OUTPUT
%      f    [optional] generated figure handle
%
%========================  CORSMED AB © 2020 ==============================
%

if nargin < 2
    mint = seq.time(1);
    maxt = seq.time(end);
elseif nargin == 2
    maxt = seq.time(end);
elseif nargin > 3
    return
end

% Expected principal time grid step (order of raster time of RF or gradients)
tres = median(seq.time(2:end)-seq.time(1:end-1)); % grid time for resampling on a full view

fig=figure;
if nargout>0
    f=fig;
end
ax=zeros(1,6);
for i=1:6
    ax(i)=subplot(3,2,i);
end
ax=ax([1 3 5 2 4 6]);   % Re-order axes as column left for RF, right for gradients
arrayfun(@(x)hold(x,'on'),ax);
arrayfun(@(x)grid(x,'on'),ax);

% Define label arrays
labels={'ADC/labels','RF mag (T)','RF/ADC ph (rad)','Gx (T/m)','Gy (T/m)','Gz (T/m)'};
arrayfun(@(x)ylabel(ax(x),labels{x}),1:6);

% populate ADC receiver position and RF info
plot(seq.time(seq.rxSignal>0), seq.rxSignal(seq.rxSignal>0), "*r", 'Parent',ax(1));
if isfield(seq,'rfScale')
    plot(seq.time, seq.rfScale, 'Parent',ax(2));
else
    plot(seq.time, abs(seq.rfmSignal), 'Parent',ax(2));
end
plot(seq.time, angle(exp(1i*angle(seq.rfmSignal)).*exp(1i*seq.rfpSignal).*exp(1i*2*pi*seq.time.*seq.rffSignal)), 'Parent',ax(3));
if isfield(seq,'swcSignal')
    % show contextually to RF magnitude trace whether any SWC were provided
    plot(seq.time(seq.swcSignal>0), seq.swcSignal(seq.swcSignal>0), "*r", 'Parent',ax(2));
    legend(ax(2),{'RF mag (T)','SWC'})
end
xlabel(ax(3),'t (s)');

% Simple quick and dirty - expand gradients from their time scale blocked into finer grid to avoid
% fake ramps due to inhomogeneous time grid adding a breakpoint
diff = seq.time(2:end)-seq.time(1:end-1);
vt = zeros(length(seq.time) + sum(diff>2*tres),1);
Gx = zeros(length(seq.time) + sum(diff>2*tres),3);
pos = 1;
brk = find(diff>2*tres);
for ii = 1 : length(seq.time)
    vt(pos) = seq.time(ii);
    Gx(pos,:) = [seq.gxSignal(ii), seq.gySignal(ii), seq.gzSignal(ii)];
    if intersect(ii,brk) > 0
        vt(pos+1) = seq.time(ii)+tres;
        pos = pos + 1;
    end
    pos = pos + 1;
end

% Plot reformatted gradient traces
plot(vt, Gx(:,1),'Parent',ax(4));
plot(vt, Gx(:,2),'Parent',ax(5));
plot(vt, Gx(:,3),'Parent',ax(6));
xlabel(ax(6),'t (s)');

% Set axis limits - link axes on same time grid (x-axis)
dispRange = [mint maxt];
arrayfun(@(x)xlim(x,dispRange),ax);
linkaxes(ax(:),'x')
h = zoom(fig);
setAxesZoomMotion(h,ax(1),'horizontal');
end
         