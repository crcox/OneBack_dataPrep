function [x] = tmp(a)
  datadir = '~/MRI/PictureOneBack/data/full';
  d = dir(datadir);
  d = d(~[d.isdir]);
  z = ~cellfun(@isempty, cellfun(@(x) regexp(x, sprintf('MRH026_4[0-9]{2}_%s_full',a)), {d.name},'Unif', 0));
  d = d(z);
  
  files = sort({d.name});
  n = numel(d);
  
  for i=1:n
    load(fullfile(datadir, files{i}),'X');
    disp(files{i})
    if i == 1
      x = generate_outlier_filters(X);
    else
      x(i,:) = generate_outlier_filters(X);
    end
  end
end

      
      