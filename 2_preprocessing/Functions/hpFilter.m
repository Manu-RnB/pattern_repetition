function lwdata = hpFilter(lwdata, subName, group)

global Cfg Log

    filerType   = Cfg.Preprocessing.HPfilter.filterType; 
    cutOff      = Cfg.Preprocessing.HPfilter.cutOff;
    order       = Cfg.Preprocessing.HPfilter.Order;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% High-Pass Butterworth filter %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    fprintf ('  -> High-Pass filter \n')

    option = struct ('filter_type', filerType,...
                     'low_cutoff', cutOff,...
                     'filter_order', order,...
                     'suffix','butt','is_save',0);
    lwdata = FLW_butterworth_filter.get_lwdata(lwdata,option);  

    fprintf ('  ==> Signal filtered \n')
    
%% Save Log file

    Log.(group).(subName).Preprocessing.HPFilter.Type       = filerType;
    Log.(group).(subName).Preprocessing.HPFilter.CuttOff    = cutOff;
    Log.(group).(subName).Preprocessing.HPFilter.Order      = order;
    
end

