clc,clear
% 30_nodc
p = mfilename('fullpath');
p = p(1:length(p)-length('run_all_op'));
cd([p,'1_30_nodc'])

    gaplist=[1,0.5,0.1,0.01]/100;
    igap= 4 ;    
    test_type='_case1_';
    diary(['op_text',test_type,num2str(igap),'_.txt'])
    %*****************************OP
    
    nj=365;
    c_reCur=0.0;
    load(['OSdata',test_type,num2str(igap)])
    
    Operation_simulation_curPunish_new_v4
    op_result(1,igap)=(c_line_lc'*Iv)/1.0e+05;
    op_result(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result(3,igap)=objsum/1.0e+05;
    op_result(4,igap)=rate*100;
    op_result(5,igap)=loadshed_energy;
    op_result(6,igap)=reCurPunish/1.0e+05;   
   
   
    save(['OSresult',test_type,num2str(igap)],'op_result')
    
    Operation_simulation_curPunish_new_v4_wet
    op_result_wet(1,igap)=(c_line_lc'*Iv)/1.0e+05;
    op_result_wet(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result_wet(3,igap)=objsum/1.0e+05;
    op_result_wet(4,igap)=rate*100;
    op_result_wet(5,igap)=loadshed_energy;
    op_result_wet(6,igap)=reCurPunish/1.0e+05;   
   
   
    save(['OSresult_wet',test_type,num2str(igap)],'op_result_wet')
    
diary off



%% 30_dc_given
cd([p,'2_30_dc'])

    test_type='_case1_';
    diary(['op_text',test_type,num2str(igap),'_2.txt'])
    
    %*****************************OP
    
%     op_result=zeros(6,length(gaplist));
    nj=365;
    c_reCur=0.0;
    load(['OSdata',test_type,num2str(igap)])
    %---transfer to OSdara_DC---
    VXdc=zeros(1,9);
    VXdc(2)=1;
    %---transfer to OSdara_DC---
    
    Operation_simulation_curPunish_DC_new_v4
    op_result(1,igap)=(c_line_lc'*Iv+c_line_lcDC'*Xv_dc)/1.0e+05;
    op_result(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result(3,igap)=objsum/1.0e+05;
    op_result(4,igap)=rate*100;
    op_result(5,igap)=loadshed_energy;
    op_result(6,igap)=reCurPunish/1.0e+05;      
   
    save(['OSresult_dc_given2',test_type,num2str(igap)],'op_result')
    
    Operation_simulation_curPunish_DC_new_v4_wet
    op_result_wet(1,igap)=(c_line_lc'*Iv+c_line_lcDC'*Xv_dc)/1.0e+05;
    op_result_wet(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result_wet(3,igap)=objsum/1.0e+05;
    op_result_wet(4,igap)=rate*100;
    op_result_wet(5,igap)=loadshed_energy;
    op_result_wet(6,igap)=reCurPunish/1.0e+05;     
   
   
    save(['OSresult_wet_dc_given2',test_type,num2str(igap)],'op_result_wet')
    
diary off

%% 50_nodc
cd([p,'3_50_nodc'])
    
    %*****************************OP
    test_type='_case3_';
%     op_result=zeros(6,length(gaplist));
    nj=365;
    c_reCur=0.0;
    load(['OSdata',test_type,num2str(igap)])
    
    Operation_simulation_curPunish_new_v4
    op_result(1,igap)=(c_line_lc'*Iv)/1.0e+05;
    op_result(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result(3,igap)=objsum/1.0e+05;
    op_result(4,igap)=rate*100;
    op_result(5,igap)=loadshed_energy;
    op_result(6,igap)=reCurPunish/1.0e+05;   
   
   
    save(['OSresult',test_type,num2str(igap)],'op_result')
    
    Operation_simulation_curPunish_new_v4_wet
    op_result_wet(1,igap)=(c_line_lc'*Iv)/1.0e+05;
    op_result_wet(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result_wet(3,igap)=objsum/1.0e+05;
    op_result_wet(4,igap)=rate*100;
    op_result_wet(5,igap)=loadshed_energy;
    op_result_wet(6,igap)=reCurPunish/1.0e+05;   
   
   
    save(['OSresult_wet',test_type,num2str(igap)],'op_result_wet')
    
diary off


%% 50_dc
cd([p,'4_50_dc'])
    
    %*****************************OP
    test_type='_case4_';
%     op_result=zeros(6,length(gaplist));
    nj=365;
    c_reCur=0.0;
    load(['OSdata',test_type,num2str(igap)])
    
    Operation_simulation_curPunish_DC_new_v4
    op_result(1,igap)=(c_line_lc'*Iv+c_line_lcDC'*Xv_dc)/1.0e+05;
    op_result(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result(3,igap)=objsum/1.0e+05;
    op_result(4,igap)=rate*100;
    op_result(5,igap)=loadshed_energy;
    op_result(6,igap)=reCurPunish/1.0e+05;      
   
    save(['OSresult',test_type,num2str(igap)],'op_result')
    
    Operation_simulation_curPunish_DC_new_v4_wet
    op_result_wet(1,igap)=(c_line_lc'*Iv+c_line_lcDC'*Xv_dc)/1.0e+05;
    op_result_wet(2,igap)=(sum(obj_mon)+sum(sum(reUpper))*c_reCur)/1.0e+05;
    op_result_wet(3,igap)=objsum/1.0e+05;
    op_result_wet(4,igap)=rate*100;
    op_result_wet(5,igap)=loadshed_energy;
    op_result_wet(6,igap)=reCurPunish/1.0e+05;     
   
   
    save(['OSresult_wet',test_type,num2str(igap)],'op_result_wet')
    
diary off

