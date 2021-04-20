%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Semi-Interactive data analysis for brain networks. 
%               For fMRI data and nodes in AAL atlas.
%                      by Isaura Oliver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Currently only analyses the networks from: MECMI, CC and PC
% 
% How it works: 
%   1. Feed the program a data set (i.e., Midnight).
%   2. Infers the Concatenate, Subject or/and Recording Networks by using CC,
%   PC and MECMI[*]. 
%   => Highly recommendended after this step to save the networks and then
%   load the inferred network to do the rest of the analysis. [Use space on
%   line XXX to load the network files in *.mat format]
%   3. Load the networks inferred in step 2.
%   4. Calculates links, degree and ranking.
%   5. Variability analysis using precision-recall, to determine how similar is the structure of the inferred networks.
%   
%
% Midnight data set: http://dx.doi.org/10.1016/j.neuron.2017.07.011
%
% [*] The MECMI method is from the paper: https://doi.org/10.1038/s41598-017-06208-w
%
% The results of this analysis are published on: https://doi.org/10.3390/e21090882
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Some Global variables
%   
% It is LABELING variable to keep track of the Midnight data set and select as nodes as are needed.
% It has the LABELS 20 nodes (DMN+FPN) in the order that is useful
% for MY Midnight dataset and not having to load the main file every time.
% It has AngL/R repeated (since both exists on DMN and FPN)
% You may ignore/change this if you are using a different atlas or more nodes.
n1 = {'MSFL'; 'MSFR'; 'ACL'; 'ACR'; 'PCL'; 'PCR'; 'AngL'; 'AngR';
       'PrecL'; 'PrecR'; 'MFL'; 'MFR'; 'IFOperL'; 'IFOperR'; 'IFTriL';
       'IFTriR'; 'IPL'; 'IPR';'AngL'; 'AngR'};
xbin = 3;
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list0 = {'Load Time-series', 'Inference Networks',  'Load Networks', ...
        'Analysis menu'};
[indx0, tf0] = listdlg( 'PromptString', 'Type Analysis', 'SelectionMode', ...
        'single', 'ListSize', [250 140], 'ListString', list0);


if indx0 ==1
        disp('Load Time-series. Useful to load if you need to check ROIs or/and do inference.');
        disp('You do not want to load every time it if you are treating big data (as my case). Use Load networks, instead.');
%       Use the path to the file (instead of uigetfile):
%       data = load('Path To your data set');
        nL = data.ROIlist; %Extracts the lists regions of interest that we need.
        %TAKE into account my structure looks like: data = dataTS and
        %ROIlist.
        
        %Give the format of the data, in this case: SubjectsxRecordingsxROIxTS. 
        %In Midnight: 10 first nodes - DMN. 10 last nodes - FPN.
        size(data.dataDMNFPN)
        %Reshape the 4-D matrix for convenience. ROIxTSxRecordingxSubject
        data0 = permute(data.dataDMNFPN, [3 4 2 1]);
        %If you need to subtract a subject, as in our case: subject 8, in the Midnight paper is marked as during the experiment suffered from drowsiness.
        %data0 = data0(:,:,:,[1:7 9:10]);
        
elseif indx0 == 2
        disp('Network inference: CC, PC, MECMI 2 and MECMI 3.');
        disp('CC uses cross correlation to infer the network, and then we threshold picked by using surrogate series.');
        disp('PC uses partial correlation to infer the network, threshold picked by surrogate series.');
        disp('MECMI uses Maximum Entropy Conditioning Mutual Information. Automatically does 2 bins and 3 bins.')

        
        list1 = {'Make surrogates', 'Inference CC networks', 'Inference PC networks', 'Inference MECMI networks'};
        [indx1, tf1] = listdlg( 'PromptString', 'Type Inference', 'SelectionMode', ...
        'single', 'ListSize', [250 140], 'ListString', list1);
        if indx1 == 1
            % Make the surroagate series/networks and find the thresholds.    
            % We use the surrogate series in order to obtain the threshold for each link in the network.
            disp('This part uses a lot memory (if big data). Once Surrogate series are done, they do the networks.')
            disp('It only saves the threshold network. If you need to save the surrogate networks of CC and PC write it output of the function "surrogateNets".');
            
            % Asks the confidence interval to threshold the surrogates and obtain the "threshold networks" for CC and PC.
            prompt = {'Enter confidence interval:'};
            dlgtitle = 'Threshold';
            definput = {'0.95'};
            dims = [1 35];
            tmp = inputdlg(prompt,dlgtitle,dims,definput);
            th = str2double(tmp{1});
            
            answer0 = questdlg('Select Surrogates needed', 'Menu', 'Concatenate', 'Subject', 'Subject-Recording', 'cancel');
            switch answer0                
                case 'Concatenate'
                    [rhoTh, prhoTh] = surrogateNets(data0,2,th); 
                    %save('pathtofile','rhoTh', 'prhoTh','-v7.3');
                case 'Subject'
                    [rhoSTh, prhoSTh] = surrogateNets(data0,3,th); 
                    %save('pathtofile','rhoSTh', 'prhoSTh','-v7.3');
                case 'Subject-Recording'
                    [rhoSRTh, prhoSRTh] = surrogateNets(data0,4,th); 
                    %save('pathtofile','rhoRSTh', 'prhoRSTh','-v7.3');
            end
                    
                
        elseif indx1 == 2
            % Asks if we need to load the threshold networks.
            answer1 = questdlg('You need to load Threshold Networks?', 'Menu', 'Yes',  'No', 'cancel');
            switch answer1
                case 'Yes' 
                    %load('path to file concatenate TH networks')
                    %load('path to file STH networks')
                    %load('path to file RSTH networks')
                case 'No'
                    disp('Okay!')
            end
            
            % First select which nodes you want to infer (for DMN, FPN or DMN+FPN)
            nodes = questdlg('Select nodes for the functional networks.', 'Menu', 'DMN',  'FPN', 'DMN+FPN', 'cancel');
            
            switch nodes
                
                case'DMN'
                    % nodes 1 to 10.
                    % We are interested on having AngL and AngR in last
                    % positions since are repated on FPN. ANd their timeseries are EXACTLY the same!
                    n = [1:6 9:10 7 8]; disp(n1(n))
                    % Time-series:
                    dataX = data0(n,:,:,:);
                case 'FPN'
                    % nodes 11 to 20.
                    % We are interested on having AngL and AngR in first positions. 
                    n = [19 20 11:18]; disp(n1(n))
                    % Time-series:
                    dataX = data0(n,:,:,:);
                    
                case 'DMN+FPN'
                    % All nodes: AngL and Ang on the middle and only once.
                    n = [1:6 9:10 7 8 11:18]; disp(n1(n))
                    dataX = data0(n,:,:,:);
            end
            
            % Infer CC networks.
            
            % concatenate network
            rhoTH = rhoTH(n, n);
            [CCc, ~] = CC_inference(dataX, rhoTH, 2); 
            % Subject networks.
            rhoSTH = rhoSTH(n, n,:);
            [CCS, ~] = CC_inference(dataX, rhoSTH, 3);
            % RS networks
            rhoRSTH = rhoRSTH(n, n,:,:);
            [CCRS, ~] = CC_inference(dataX, rhoRSTH, 4);
                        
            switch nodes
                case 'DMN'
                    CCc_D = CCc; CCS_D = CCS; CCRS_D = CCRS;
                    
                case 'FPN'
                    CCc_F = CCc; CCS_F = CCS; CCRS_F = CCRS;
                case 'DMN+FPN'
                    CCc_DF = CCc; CCS_DF = CCS; CCRS_DF = CCRS;
            end
            
            answer1a = questdlg('Save the networks or continue to infer the other functional networks?', 'Menu', 'Save',  'Continue', 'cancel');
            switch answer1a
                case 'Save' 
                    %If you run the inference routine for the three networks you can save them all on one file:
                    %save('pathtofile','CCc_D', 'CCc_F', 'CCc_DF', 'CCS_D', 'CCS_F', 'CCS_DF', 'CCRS_D', 'CCRS_F', 'CCRS_DF','-v7.3');
                case 'Continue'
                    disp('Then rerun the program and select the next functional network you want!')
            end
            
            
        elseif indx1 == 3
            % Inference of PC networks.
            % concatenate network
            [PCc, ~] = PC_inference(dataX, prhoTH, 2); 
            % Subject networks.
            [PCS, ~] = PC_inference(dataX, prhoSTH, 3);
            % RS networks.
            [PCRS, ~] = PC_inference(dataX, prhoRSTH, 4);
            
            switch nodes
                case 'DMN'
                    PCc_D = PCc; PCS_D = PCS; PCRS_D = PCRS;
                    
                case 'FPN'
                    PCc_F = PCc; PCS_F = PCS; PCRS_F = PCRS;
                case 'DMN+FPN'
                    PCc_DF = PCc; PCS_DF = PCS; PCRS_DF = PCRS;
            end
            
            answer1aa = questdlg('Save the networks or continue to infer the other functional networks?', 'Menu', 'Save',  'Continue', 'cancel');
            switch answer1aa
                case 'Save' 
                    %If you run the inference routine for the three networks you can save them all on one file:
                    %save('pathtofile','PCc_D', 'PCc_F', 'PCc_DF', 'PCS_D', 'PCS_F', 'PCS_DF', 'PCRS_D', 'PCRS_F', 'PCRS_DF','-v7.3');
                case 'Continue'
                    disp('Then rerun the program and select the next functional network you want!')
            end
        elseif indx1 == 4
            % MECMI inference in 3 bins, since we want to take into account
            % the non-linearity of our data.
            L = size(dataX,1);
            N = 10; % it will make an ensamble of 10 networks which we will pick the one with the max(MaxEntropy).
            
            % concatenate network
            [M3, Sm] = MECMI_inference(dataX, L, N, xbin, 2); 
            % Subject networks.
            [M3S, SmS] = MECMI_inference(dataX, L, N, xbin, 3);
            % RS networks.
            [M3RS, SmRS] = MECMI_inference(dataX, L, N, xbin, 4);
            
            switch nodes
                case 'DMN'
                    M3_D = M3; M3S_D = M3S; M3RS_D = M3RS;
                    
                case 'FPN'
                    M3_F = M3; M3S_F = M3S; M3RS_F = M3RS;
                case 'DMN+FPN'
                    M3_DF = M3; M3S_DF = M3S; M3RS_DF = M3RS;
            end
        end
        
        % save the networks as:
        %save('Path to the folder and filename (it saves with *.mat format).','namesofALLnetworks', '-v7.3');
elseif indx0 == 3
        disp('Load the functional networks for further analysis.');
        disp('For this analysis we require 3 types of Networks adquired in 3 different similitude measures and 3 functional networks.');
        disp('Functional networks: DMN, FPN and DMN+FPN.')
        disp('Types of inference: CC, PC, MECMI');
        disp('Types of networks: concatenated networks, Subject networks and Recording-Subject networks');
        % To load the Networks, since they may be a few files and manual
        % search is tiring use the load('pathToFile'):
        
        %Concatenate Networks.
        %data1 = load('Load PC Networks for DMN, FPN, DMN+FPN.');
        %data2 = load('Load PC Networks for DMN, FPN, DMN+FPN.');
        %data3 = load('Load MECMI in 3 bins Networks for DMN, FPN, DMN+FPN.');
                %if you also need to load MECMI in 2 bins add the
                %correspondent file and remember to also write the networks
                %below!
        
        xbin = 3; %Tells us that we are using MECMI 3.        
        % Tell which Functional networks you want to load
        answer1a = questdlg('Select nodes you want', 'Menu', 'DMN', 'FPN', 'All nodes', 'cancel');
        switch answer1a
        
            case 'DMN'
                CC = data1.CC_D; PC = data2.PC_D; Mc = data3.M3_D;
                CCS = data1.CCS_D; PCS = data2.PCS_D; MS = data3.MS3_D;
                CCRS = data1.CCRS_D; PCRS = data2.PCRS_D; MRS = data3.MSR3_D;    
                    
            case 'FPN'
                CC = data1.cc_F; PC = data2.PC_F; Mc = data3.M3_F;
                CCS = data1.CCS_F; PCS = data2.PCS_F; MS = data3.MS3_F;
                CCRS = data1.CCRS_F; PCRS = data2.PCRS_F; MRS = data3.MSR3_F;

            case 'All nodes'
                CC = data1.CC_DF; PC = data2.PC_DF; Mc = data3.M3_DF;
                CCS = data1.CCS_DF; PCS = data2.PCS_DF; MS = data3.MS3_DF;
                CCRS = data1.CCRS_DF; PCRS = data2.PCRS_DF; MRS = data3.MSR3_DF;
        end
    

elseif indx0 ==4
%   'Analysis menu'    
    list = {'Degree & Link distribution', 'Ranking',  'Variability', ...
        'Variability vs MECMI3', 'Ranking Lists', 'Intra-links variability', 'Inter-links variability', ...
        'Indirect links', 'Labeling', 'Empty'};
    [indx, tf] = listdlg( 'PromptString', 'Type Analysis', 'SelectionMode', ...
        'single', 'ListSize', [250 140], 'ListString', list);
    mats = questdlg('Matrices you want to use?', 'Menu', 'MECMI', 'PC', 'CC', 'cancel');
        switch mats
            case 'MECMI'
                A1 = Mc;
                A2 = MS;
                A3 = MRS;
            case 'PC'
                A1 = PC; 
                A2 = PCS; 
                A3 = PCRS;
            case 'CC'
                A1 = CC;
                A2 = CCS;
                A3 = CCRS;     
        end
        
    if indx == 1
    %   1. Degree and link distribution
        [degC, degS, xT, yT] = Degree(A1,A2);
        disp(degC)
        disp(degS)
        % Keep the degrees
        switch mats
            case 'MECMI'
                degM = degC;
                degMS = degS;
                if xbin == 3
                    degM3 = degM;
                    degMS3 = degMS; 
                end
            case 'PC'
                degP = degC;
                degPS = degS;
            case 'CC'
                degCC = degC;
                degCCS = degS;
        end
        
        switch answer1a
        	case 'DMN'
                switch mats
                	case 'MECMI'
                        X = 'M3D';
                    case 'PC'
                        X = 'PCD';
                    case 'CC'
                        X = 'CCD';
                end
            case 'FPN'
                switch mats
                    case 'MECMI'
                    	X = 'M3F';
                    case 'PC'
                        X = 'PCF';
                    case 'CC'
                        X = 'CCF';
                end
            case 'DNM+FPN'
                switch mats
                    case 'MECMI'
                        X = 'M3DF';
                    case 'PC'
                        X = 'PCDF';
                    case 'CC'
                        X = 'CCDF';
                end
        end
        %If you want to save the histogram of the Subject networks.
        %prompt = {'Folder','Distribution for Subject networks.'};
        %title = 'Save the files';
        %numlines = 1;
        %defaultname = {['degreeS_' X '.dat']};
        %name = inputdlg(prompt,title,numlines,defaultname);
        %f = fullfile('Path to foldersave', 'degree', name{1});
        %dlmwrite(f,[xThis' yThis'], 'delimiter', ' ');
        
        
    elseif indx == 2
        
        switch answer1a
            case 'DMN'
                switch mats
                	case 'MECMI'
                    	X = 'M3D';
                    case 'PC'
                    	X = 'PCD';
                    case 'CC'
                    	X = 'CCD';
                end
            case 'FPN'
                switch mats
                	case 'MECMI'
                    	X = 'M3F';
                    case 'PC'
                    	X = 'PCF';
                    case 'CC'
                        X = 'CCF';
                end
            case 'DMN+FPN'
                switch mats
                	case 'MECMI'
                    	X = 'M3DF';
                    case 'PC'
                    	X = 'PCDF';
                    case 'CC'
                    	X = 'CCDF';
                end
        end
        
    %   2. Ranking        
        for i = 1:size(A2,3)
            % concatenate vs subjects
            ranked = Ranking(A1,A2(:,:,i));
            disp('To save the files uncomment.')
            % save the files
            %f = fullfile('Path', 'ToFolder', 'degree', ['rank_S' num2str(i) '_' X '.dat']);
            %dlmwrite(f,ranked, 'delimiter', ' ');
            
                % subject vs subject
                for j = 1:size(A2,3)
                    rankS = Ranking(A2(:,:,i), A2(:,:,j));
                %   Save the files    
                %   f = fullfile('Path', 'ToFolder', 'degree', ['rank_S_' num2str(i) 'vsS_' num2str(j) '_' X '.dat']);
                %   dlmwrite(f,ranked, 'delimiter', ' ');
                end
            
                % subject network X vs its own recordings network
                for j = 1:size(A3,3)
                    rankRS = Ranking(A2(:,:,i), A3(:,:,j,i));
                %   Save the files        
                %   f = fullfile('Path', 'ToFolder', 'degree', ['rank_S_' num2str(i) 'vsR_' num2str(j) '_' X '.dat']);
                %   dlmwrite(f,ranked, 'delimiter', ' ');
                end
        
        end
%       2.b. Ranks Method vs Method: concatenate vs concatenate
        matis = questdlg('Select comparision? For Ranking ', 'Menu', 'MECMI vs PC', ...
            'MECMI vs CC', 'PC vs CC', 'cancel');
        switch matis
            case 'MECMI vs PC'
                AM = Mc;
                AX = PC; 
                % M3 concatenate vs PC concatenate
                rankedvs3 = Ranking(AM,AX);
                switch answer1a
                    case 'DMN'
                        Y = 'M3D';
                        X = 'PCD';
                    case 'FPN'
                        Y = 'M3F';
                        X = 'PCF';                          
                    case 'DMN+FPN'
                        Y = 'M3DF';
                        X = 'PCDF';
                end     
            
            %save the files.
            prompt = {'folder', 'Ranking for M3 vs PC'}; 
            title = 'Name of files for saving';
            numlines = 1;
            defaultname = {['ranking_' Y '_vs_' X '.dat']};
            name = inputdlg(prompt,title,numlines,defaultname);
            %f = fullfile('Path', 'ToFolder', 'degreeRank', name{1});
            %dlmwrite(f,rankedvs3, 'delimiter', ' ');

            case 'MECMI vs CC'
                AM = Mc;
                AX = CC;
                
                rankedvs3 = Ranking(AM,AX);
                
                switch answer1a
                    case 'DMN'
                        Y = 'M3D';
                        X = 'CCD';
                    case 'FPN'
                        Y = 'M3F';
                        X = 'CCF';                          
                    case 'DMN+FPN'
                        Y = 'M3DF';
                        X = 'CCDF';
                end     
            
                %save the files.
                prompt = {'folder', 'Ranking for M3 vs CC'}; 
                title = 'Name of files for saving';
                numlines = 1;
                defaultname = {['ranking_' Y '_vs_' X '.dat']}; 
                name = inputdlg(prompt,title,numlines,defaultname);
                %f = fullfile('Path', 'ToFolder', 'degreeRank', name{1});
                %dlmwrite(f,rankedvs3, 'delimiter', ' ');
                
                          
            case 'PC vs CC'
                AM = PC; %change to PC if want to do bins or P99 for surrogates
                AX = CC;
                
                rankedvs = Ranking(AM,AX);
                
                switch answer1a
                    case 'DMN'
                        Y = 'PCD';
                        X = 'CCD';
                    case 'FPN'
                        Y = 'PCF';
                        X = 'CCF';                          
                    case 'DMN+FPN'
                        Y = 'PCDF';
                        X = 'CCDF';
                end 
                
                %save the files 
                prompt = {'folder','Ranking for PC vs CC'};
                title = 'Name of files for saving';
                numlines = 1;
                defaultname = {['ranking_' Y '_vs_' X '.dat']};
                name = inputdlg(prompt,title,numlines,defaultname);
                %f = fullfile('Path', 'ToFolder', 'degreeRank', name{1});
                %dlmwrite(f,rankedvs, 'delimiter', ' ');
                
        end
        
        
    elseif indx == 3
    %   3. Variability on same network
        [degC, degS, ~, ~] = Degree(A1,A2);
        
        typix = questdlg('Variability for full matrix or some links only?', 'Menu', ...
            'Full matrix', 'Some links', 'cancel');
        typix2 = questdlg('Variability for full matrix or some links only?', 'Menu', ...
             'Concatenate vs Subject', 'Subject vs Subject', 'Subj vs Record', 'cancel');
        switch typix
            case 'Full matrix'
            % 1. Full matrix
                prot = 'Enter title';
                title = 'Enter title';
                dims = [1 35];
                definpt = {'FPN'};
                titl = inputdlg(prot, title, dims, definpt);
                % saves ?
                switch answer1a
                    case 'DMN'
                        switch mats
                            case 'MECMI'
                                X = 'M3D';
                            case 'PC'
                                X = 'PCD';
                            case 'CC'
                                X = 'CCD';
                        end
                    case 'FPN'
                        switch mats
                            case 'MECMI'
                                X = 'M3F';
                            case 'PC'
                                X = 'PCF';
                            case 'CC'
                                X = 'CCF';
                        end                         
                    case 'DMN+FPN'
                        switch mats
                            case 'MECMI'
                                X = 'M3DF';
                            case 'PC'
                                X = 'PCDF';
                            case 'CC'
                                X = 'CCDF';
                        end
                end     
                
             switch typix2   
               %concatenate vs subjects
             case 'Concatenate vs Subject'
                matVar = Variability(A1,A2, degC, titl);  
                %f = fullfile('Path', 'ToFolder', 'variability', ['var_Full_' X '.dat']);
                %dlmwrite(f,matVar, 'delimiter', ' ');
             case 'Subject vs Subject'    
                %subjects vs subjects
                for i = 1:size(A2,3)
                    matVar = Variability(A2(:,:,i),A2, degS(i), titl);
                    %f = fullfile('Path', 'ToFolder', 'variability', ['var_Full_S' num2str(i) '_' X '.dat']);
                    %dlmwrite(f,matVar, 'delimiter', ' ');
                end
             case 'Subj vs Record'    
                %subjects vs recordings
                for i = 1:size(A2,3)
                    % subject X vs all his recordings.
                    matVar = Variability(A2(:,:,i),A3(:,:,:,i), degS(i), titl);
               
                    % saves ?                   
                    %f = fullfile('Path', 'ToFolder', variability', ['var_S' num2str(i) 'vsRX_' X '.dat']);
                    %dlmwrite(f,matVar, 'delimiter', ' ');
                end
             end                   
                
            case 'Some links'
                % 2. Only some links
                promt = {'Enter number of links', 'Enter title graph'};
                titlo = 'Input';
                dims = [1 25];
                definput = {'5', 'FPN'};
                xarg = inputdlg(promt, titlo, dims, definput);
                x = str2double(xarg{1});
                
                switch answer1a
                    case 'DMN'
                        switch mats
                            case 'MECMI'
                                X = 'M3D';
                            case 'PC'
                                X = 'PCD';
                            case 'CC'
                                X = 'CCD';
                        end
                    case 'FPN'
                        switch mats
                            case 'MECMI'
                                X = 'M3F';
                            case 'PC'
                                X = 'PCF';
                            case 'CC'
                                X = 'CCF';
                        end                         
                    case 'DMN+FPN'
                        switch mats
                            case 'MECMI'
                                X = 'M3DF';
                            case 'PC'
                                X = 'PCDF';
                            case 'CC'
                                X = 'CCDF';
                        end
                end 
                deg = str2double(xarg{1});
               
                % concatenate
                tmp = triu(A1,1); 
                [srt, lis] = sort(tmp(:),'descend');
                A1Xb = A1.*(A1 >= srt(x));
            
                % subjects
                A2Xb = zeros(size(A2,1), size(A2,2), size(A2,3)); 
                for k = 1:size(A2,3)
                    tmp = triu(A2(:,:,k),1); 
                    srtp = sort(tmp(:),'descend'); 
                    A2Xb(:,:,k) = A2(:,:,k).*(A2(:,:,k) >= srtp(x));
                end
                
                % recordings
                A3Xb = zeros(size(A3,1), size(A3,2), size(A3,3), size(A3,4)); 
                for k = 1:size(A3,3)
                    for l = 1:size(A3,4)
                        tmp = triu(A3(:,:,k,l),1); 
                        srtp = sort(tmp(:),'descend'); 
                        A3Xb(:,:,k,l) = A3(:,:,k,l).*(A3(:,:,k,l) >= srtp(x));
                    end
                end
                
               switch typix2   
                    case 'Concatenate vs Subject'
                        matVar = Variability(A1Xb,A2Xb, deg, xarg{2});
                        %f = fullfile('PathToFolder', 'variability', ['links' xarg{1} '/var_' xarg{1} '_' X '.dat']);
                        %dlmwrite(f,matVar, 'delimiter', ' ');
                        
                    case 'Subject vs Subject'    
                        %subjects vs subjects                 
                        for i = 1:size(A2,3)
                            matVar = Variability(A2Xb(:,:,i),A2Xb, deg, titl);
                            %f = fullfile('PathToFolder', 'variability', ['links' xarg{1} '/var_' xarg{1} '_S' num2str(i) '_' X '.dat']);
                            %dlmwrite(f,matVar, 'delimiter', ' ');
                        end
                    case 'Subj vs Record'    
                    %subjects vs recordings
                    for i = 1:size(A2,3)
                        matVar = Variability(A2Xb(:,:,i),A3Xb(:,:,:,i), deg, titl); 
                                        
                        % saves ?                   
                        %f = fullfile('PathToFolder', variability', ['links' xarg{1} '/var_' xarg{1} '_S' num2str(i) 'vsRX_' X '.dat']);
                        %dlmwrite(f,matVar, 'delimiter', ' ');
                    end
               end
               
        end
        
        elseif indx == 4
    %   4. Variability study comparing CC and PC to MECMI 3
        [degC, degS, ~, ~] = Degree(A1,A2);
        switch answer1a
        	case 'DMN'
                switch mats
                	case 'PC'
                        X = 'PCD';
                    case 'CC'
                        X = 'CCD';
                end
            case 'FPN'
                switch mats
                	case 'PC'
                    	X = 'PCF';
                    case 'CC'
                        X = 'CCF';
                end
            case 'DMN+FPN'
                switch mats
                	case 'PC'
                    	X = 'PCDF';
                    case 'CC'
                    	X = 'CCDF';
                end   
        end
        
        
        typix = questdlg('Variability for full matrix or some links only?', 'Menu', ...
            'Full matrix', 'Some links', 'cancel');
        switch typix
                  
            case 'Full matrix'
            % 1. Full matrix
                prot = 'Enter title';
                title = 'Enter title';
                dims = [1 35];
                definpt = {'FPN'};
                titl = inputdlg(prot, title, dims, definpt);
                
                MC3 = Mc;
                [degM3, degS, ~, ~] = Degree(MC3,A2);
                matVar = Variability(MC3,A2, degM3, titl);
                
                % save
                prompt = {'folder','Variability for X vs MECMI3 Full'};
                title = 'Name of files for saving';
                numlines = 1;
                defaultname = {['var_FULL_MECMI3_' X '.dat']};
                name = inputdlg(prompt,title,numlines,defaultname);
                %f = fullfile('PathToFolder', 'variability', name{2});
                %dlmwrite(f,matVar, 'delimiter', ' ');
            
                
            case 'Some links'
            % 2. Only some nodes
            % number of nodes to do analysis
                promt = {'Enter number of links', 'Enter title graph'};
                titlo = 'Input';
                dims = [1 25];
                definput = {'5', 'FPN'};
                xarg = inputdlg(promt, titlo, dims, definput);
                
                x = str2double(xarg{1});
                MC3 = Mc; %only when loading matrices
                
                % concatenate
                tmp = triu(MC3,1); 
                [srt, lis] = sort(tmp(:),'descend');
                A1X = MC3.*(MC3 >= srt(x));
            
                % subjects
                A2X = zeros(size(A2,1), size(A2,2), size(A2,3)); 
                for k = 1:size(A2,3)
                    tmp = triu(A2(:,:,k),1); 
                    srtp = sort(tmp(:),'descend'); 
                    A2X(:,:,k) = A2(:,:,k).*(A2(:,:,k) >= srtp(x));
                end
                degM3Sp = str2double(xarg{1});
                matVar = Variability(A1X,A2X, degM3Sp, xarg{2});
                
                % subjects and recordings
                
                % saves ?               
                prompt = {'folder', 'Variability for X vs MECMI3 Partial'};
                title = 'Name of files for saving';
                numlines = 1;
                defaultname = {['var_' xarg{1} '_MECMI3_' X '.dat']};
                name = inputdlg(prompt,title,numlines,defaultname);
                %f = fullfile('PathToFolder', 'variability/links5', name{2});
                %dlmwrite(f,matVar, 'delimiter', ' ');
        end
        elseif indx == 5
        % Get the strongest links
        switch answer1a
            case 'DMN'
                A1ND = 0; A2ND = zeros(size(A2));
            case 'FPN'
                A1ND = 0; A2ND = zeros(size(A2));
            case 'DMN+FPN'
                A1ND = NonBlockDiag(A1); 
                A2ND = zeros(size(A2,1)-2, size(A2,2)-2, size(A2,3));
                for i = 1:size(A2,3)   
                    A2ND(:,:,i) = NonBlockDiag(A2(:,:,i));
                end         
        end
            
        % concatenate vs subject
        ranked = zeros(5,5,size(A2,3));
        indxS5 = zeros(5,size(A2,3));
        for i = 1:size(A2,3)
            [ranked(:,:,i), indxS5(:,i)] = RankingSub(A1,A2(:,:,i),A1ND,A2ND(:,:,i),5);
        end
        % subject vs subject
        ranksS = zeros(5, 5, size(A2,3));
        indxS5S = zeros(5,size(A2,3));
        for i = 1:size(A2,3)
            for j = 1:size(A2,3)
                [ranksS(:,:,j,i), indxS5S(:,i)] = RankingSub(A2(:,:,i),A2(:,:,j),A2ND(:,:,i),A2ND(:,:,j),5);
            end
        end
            
        % now tell me the index found from which two nodes correspond:
        tableRCPS = zeros(size(indxS5,1), 3, size(A2,3));
        for i = 1:size(A2,3) 
            [r, c, rankPos] = locNodes(indxS5S(:,i), A2(:,:,i));
            mate = [rankPos r c];
            tableRCPS(:,:,i) = sortrows(mate,1);
        end
            
        Roi = n1([1:6 9 10 7 8 11:18]);
        RoiAl = n1([1:6 9 10 7 8 11:18]); 
        t1 = cell(size(tableRCPS,1), 2, size(tableRCPS,3));
        t2 = cell(size(tableRCPS,1), 2, size(tableRCPS,3));
        for k = 1:size(tableRCPS,1)
            for j = 2:size(tableRCPS,2)
                for i = 1:size(tableRCPS,3)
                    t1(k,j,i) = cellstr(Roi{tableRCPS(k,j,i)});
                    t2(k,j,i) = cellstr(RoiAl{tableRCPS(k,j,i)});     
                end
            end
        end          
        %Make a table with it!   
        t1 = cell(size(tableRCPS,1), 2, size(tableRCPS,3));
        t2 = cell(size(tableRCPS,1), 2, size(tableRCPS,3));
        for k = 1:size(tableRCPS,1)
            for j = 2:size(tableRCPS,2)
                for i = 1:size(tableRCPS,3)
                    t1(k,j,i) = cellstr(Roi{tableRCPS(k,j,i)});
                    t2(k,j,i) = cellstr(RoiAl{tableRCPS(k,j,i)});
                end
            end
        end
            
            
        % subject vs recording (not useful but more lists if you like)
        ranks = zeros(5, 2, size(A3,3), size(A3,4));
        indxS5 = zeros(5,size(A3,4));
        for i = 1:size(A3,4)
            for j = 1:size(A3,3)
                [ranks(:,:,j,i), indxS5(:,i)] = RankingSub(A2(:,:,i),A3(:,:,j,i),5);
            end
        end
           
        % save            
        switch answer1a
            case 'DMN'
              switch mats
              	case 'MECMI'
                	X = 'M3D';
                    tM3 = t2;
                case 'PC'
                	X = 'PCD';
                    tPC = t2;
                case 'CC'
                	X = 'CCD';
                    tCC = t2;
              end
                
              case 'FPN'
                switch mats
                    case 'MECMI'
                        X = 'M3F';
                        tM3 = t2;
                    case 'PC'
                        X = 'PCF';
                        tPC = t2;
                    case 'CC'
                        X = 'CCF';
                         tCC = t2;
                end
               case 'DMN+FPN'
                   switch mats
                    case 'MECMI'
                    	X = 'M3DF';
                        tM3 = t2;
                    case 'PC'
                        X = 'PCDF';
                        tPC = t2;
                    case 'CC'
                        X = 'CCDF';
                        tCC = t2;
                   end
        end
        %save(['PathToFolder\rankLists\table_CvsS_' X '2.mat'], 't1');
        %save(['PathToFolder\rankLists\table_Alias_CvsS_' X '2.mat'], 't2');       
        
    elseif indx == 6       
        % 1. Non-block diagonal.
        switch answer1a
            case 'DMN+DMN'
                A1ND = NonBlockDiag(A1);
                A2ND = zeros(size(A2,1)-2,size(A2,2)-2, size(A2,3));
                for i = 1:size(A2,3)
                    A2ND(:,:,i) = NonBlockDiag(A2(:,:,i));
                end 
                A3ND = zeros(size(A3,1)-2, size(A3,2)-2, size(A3,3), size(A3,4));
                for i = 1:size(A3,3)
                    for j = 1:size(A3,4)
                        A3ND(:,:,i,j) = NonBlockDiag(A3(:,:,i,j));
                    end
                end
                
                switch mats
                    case 'MECMI'
                        X = 'M3';
                    case 'PC'
                        X = 'PC';
                        PCND = A1ND;
                    case 'CC'
                        X = 'CC';
                        CCND = A1ND;
                end
        end
        % get the degree
        [deg, degS, ~, ~] = Degree(A1ND,A2ND);
        
        % 2. Variability intra-links on the concatenates
        prot = 'Enter title';
        title = 'Enter title';
        dims = [1 35];
        definpt = {'DMN+FPN'};
        titl = inputdlg(prot, title, dims, definpt);
        matVar = Variability(A1ND,A2ND, deg, titl);
                
        prompt = {'folder','Intra-links variability'};
        title = 'Name of files for saving';
        numlines = 1;
        defaultname = {['var_FULL_' X '.dat']};
        name = inputdlg(prompt,title,numlines,defaultname);
        %f = fullfile('PathToFolder', 'intLinks', name{2});
        %dlmwrite(f,matVar, 'delimiter', ' ');
        clear matVar
        
        %Variability intra-links on Subject networks
        prot = 'Enter title';
        title = 'Enter title';
        dims = [1 35];
        definpt = {'DMN+FPN'};
        titl = inputdlg(prot, title, dims, definpt);
                
        pres = zeros(size(A2ND,3), size(A2ND,3)); rec = zeros(size(A2ND,3), size(A2ND,3));
        for i = 1:size(A2ND,3)
            [matVar, pres(i,:), rec(i,:)] = Variability(A2ND(:,:,i), A2ND, degS(i), titl);
            %f = fullfile('PathToFolder', ['var_S' num2str(i) '_' X '.dat']);
            %dlmwrite(f,matVar, 'delimiter', ' ');
        end
        
        %clear matVar
        %Variability intra-links of Subject network against their own
        %recordings.
        prot = 'Enter title';
        title = 'Enter title';
        dims = [1 35];
        definpt = {'DMN+FPN'};
        titl = inputdlg(prot, title, dims, definpt);
                
        for i = 1:size(A2ND,3)
            [matVar,~,~] = Variability(A2ND(:,:,i),A3ND(:,:,:,i), degS(i), titl);
            %f = fullfile('PathToFolder', ['var_S' num2str(i) 'vsRX_' X '.dat']);
            %dlmwrite(f,matVar, 'delimiter', ' ');
        end
        
	elseif indx == 7
        % 1. Get the block diagonal
        switch answer1a
        	case 'DMN+DMN'
                A1ND = BlockDiag(A1);
                A2ND = zeros(size(A2,1)-2,size(A2,2)-2, size(A2,3));
                for i = 1:size(A2,3)
                    A2ND(:,:,i) = BlockDiag(A2(:,:,i));
                end 
                A3ND = zeros(size(A3,1)-2,size(A3,2)-2, size(A3,3), size(A3,4));
                for i = 1:size(A3,3)
                    for j = 1:size(A3,4)
                        A3ND(:,:,i,j) = NonBlockDiag(A3(:,:,i,j));
                    end 
                end 
                
                switch mats
                    case 'MECMI'
                        X = 'M3';
                    case 'PC'
                        X = 'PC';
                    case 'CC'
                        X = 'CC';
                end
                
        end
            % degree
            [deg, ~, ~, degS, ~, ~, ~, ~] = Degree(A1ND,A2ND);
            % 2. do the variability against the concatenate
            prot = 'Enter title';
            title = 'Enter title';
            dims = [1 35];
            definpt = {'DMN+FPN'};
            titl = inputdlg(prot, title, dims, definpt);
            matVar = Variability(A1ND,A2ND, deg, titl);
                
            prompt = {'folder','Variability for inter-links.'};
            title = 'Name of files for saving';
            numlines = 1;
            defaultname = {['var_full_' X '.dat']};
            name = inputdlg(prompt,title,numlines,defaultname);
            %f = fullfile('PathToFolder', 'intLinks/block', name{1});
            %dlmwrite(f,matVar, 'delimiter', ' ');
            %clear matVar
            
            % Inter-link variability Subject vs Subject
            prot = 'Enter title';
            title = 'Enter title';
            dims = [1 35];
            definpt = {'FPN'};
            titl = inputdlg(prot, title, dims, definpt);
                
            pres = zeros(size(A2ND,3), size(A2ND,3)); rec = zeros(size(A2ND,3), size(A2ND,3));
            for i = 1:size(A2ND,3)
                [matVar, pres(i,:), rec(i,:)] = Variability(A2ND(:,:,i), A2ND, degS(i), titl);
                %f = fullfile('PathToFolder',  'intLinks/block', ['var_S' num2str(i) '_' X '.dat']);
                %dlmwrite(f,matVar, 'delimiter', ' ');
            end
            
            %clear matVar
            
            %Variability inter-links Subject vs its own recordings.
            prot = 'Enter title';
            title = 'Enter title';
            dims = [1 35];
            definpt = {'FPN'};
            titl = inputdlg(prot, title, dims, definpt);
                
            for i = 1:size(A2,3)
                matVar = Variability(A2(:,:,i),A3(:,:,:,i), degS(i), titl); 
     
                %f = fullfile('PathTofolder', 'intLinks/block', ['var_S' num2str(i) 'vsRX_' X '.dat']);
                %dlmwrite(f,matVar, 'delimiter', ' ');
                
            end
            
    elseif indx == 8
        
        %Determines if the triangles found in the CC DMN+FPN concatenate network are
        %direct or indirect by using PC concatenate network. 
        disp('This analysis needs to be run after you have got the intra-links networks.')
        %Uses the concatenate with 18 nodes networks (DMN+FPN) and intra-links 16 nodes networks.
        A1 = PC;
        A = CC;
        
        [Y, bL, xcor, ycor, ~] = Ranking_v2(A1,A,PCND,CCND);
        
        %use xcor and ycor to find the nodes of the indirect links.
        nL = n1([1:6 9:18]);
        nxcor = nL(xcor); nycor = nL(ycor); 
        %inLi = [xcor, ycor];
        
    elseif indx == 9
        % 5 strongest links in the concatenate and subject networks.
        disp('Uses files table_Alias_CvsS_X.dat to label the five strongest links in the networks: DMN, FPN and DMN+FPN.')        
        
        %tCC = load('File table_Alias_CvsS_CCD.dat)
        %tPC = load('File table_Alias_CvsS_PCD.dat)
        %tM3 = load('File table_Alias_CvsS_M3D.dat)
        [ClD, C1sorD] = linkFreqLabel(tCC,tPC,tM3);
        
        %tCC = load('File table_Alias_CvsS_CCF.dat)
        %tPC = load('File table_Alias_CvsS_PCF.dat)
        %tM3 = load('File table_Alias_CvsS_M3F.dat)
        [C1F, C1sorF] = linkFreqLabel(tCC,tPC,tM3);
        
        %tCC = load('File table_Alias_CvsS_CCDF.dat)
        %tPC = load('File table_Alias_CvsS_PCDF.dat)
        %tM3 = load('File table_Alias_CvsS_M3DF.dat)
        [C1, C1sor] = linkFreqLabel(tCC,tPC,tM3);
        
        
    
    elseif indx == 10
        
        
        
        
    end
        
end
