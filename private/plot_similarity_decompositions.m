function plot_similarity_decompositions(embedding)
    figure(1);
    subplot(1,5,1)
    edist = squareform(pdist(embedding,'euclidean'));
    plot(eig(edist))
    xlim([0,20])
    title('Euclidean Distance');
    
    subplot(1,5,2)
    cdist = squareform(pdist(embedding,'cosine'));
    plot(eig(cdist))
    xlim([0,20])
    title('Cosine Distance');
    
    subplot(1,5,3)
    pcorr = corr(embedding','type','Pearson');
    plot(real(eig(pcorr)))
    xlim([0,20])
    title('Pearson Correlation');
    
    subplot(1,5,4)
    scorr = corr(embedding','type','Spearman');
    plot(real(eig(scorr)))
    xlim([0,20])
    title('Spearman Correlation');
    
    subplot(1,5,5)
    covmat = cov(embedding');
    plot(real(eig(covmat)))
    xlim([0,20])
    title('Covariance');
end