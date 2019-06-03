function spikes = get_tf(cells, cid, trial)
% extract the required trial of a cell from temporal frequency
%     input:
%        cells = the data set
%        cid = the required cell ID
%        trial = the required trial
%     output:
%        trial = a 1*n vector of spikes

spikes = squeeze(cells(cid).tf_spikes(:,2,:))';
spikes = spikes(:);
spikes = spikes{trial};

end