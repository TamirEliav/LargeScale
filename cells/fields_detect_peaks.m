function fields = fields_detect_peaks(FR_map,prm)

warning('off', 'signal:findpeaks:largeMinPeakHeight');

PSTH_Nan2zero = FR_map.PSTH; PSTH_Nan2zero(isnan(PSTH_Nan2zero )) = 0; %changing nans to zeros so that border peaks could be found. 
[pks,locs_IX] = findpeaks(PSTH_Nan2zero,'minPeakHeight',prm.fields.FR_thr);
locs = FR_map.bin_centers(locs_IX);
fields = repelem(struct,length(locs_IX));
[fields(:).loc_IX] = disperse(locs_IX);
[fields(:).loc]    = disperse(locs);
[fields(:).peak]   = disperse(pks);

warning('on', 'signal:findpeaks:largeMinPeakHeight');

end