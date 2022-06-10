function [] = setup_miji()
    % @author: pdzialecka
    
    %% Add FIJI / MIJ paths
    % ADJUST these depending on local directories
    javaaddpath 'C:\Program Files\MATLAB\R2021a\java\jar\mij.jar'
    javaaddpath 'C:\Users\Pat\Desktop\Fiji.app\jars\ij-1.53q.jar'
    addpath(genpath('C:\Users\Pat\Desktop\Fiji.app\scripts'));
    
end
