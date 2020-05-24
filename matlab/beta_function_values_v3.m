datasets = {'/homes/gullo/ProbabilisticCores/beta_function_values/Flickr.txt.graph.part.3_0_beta_matlab.txt'; '/homes/gullo/ProbabilisticCores/beta_function_values/DBLP_mu=5.txt_beta_matlab.txt'; '/homes/gullo/ProbabilisticCores/beta_function_values/biomine_probgraph.txt_beta_matlab.txt'; '/homes/gullo/ProbabilisticCores/beta_function_values/biomine_probgraph.txt.graph.part.20_3_beta_matlab.txt'};
mineta = 0.00001;

for i=1:size(datasets),
	dat = char(datasets(i))
	fid = fopen(dat);
	A = fscanf(fid,'%g %g %g',[3 inf]);
	A = A';
	fclose(fid);
	
	datout = [dat '_output'];
	fid = fopen(datout,'w');
	
	for i=1:size(A,1),
		v = A(i,1);
		p = A(i,2);
		dv = A(i,3);
		betas = [v p dv];
		betas2 = [v p dv];
		s = '%d\t%g\t%d';
		k = 0;
		stop = false;
		while (k<=dv) && (stop==false)
			x = p;
			z = k;
			w = dv-k+1;
			w2 = dv-k;
			I = betainc(x,z,w);
			if I < mineta
				stop = true;
			end
			I2 = betainc(x,z,w2);
			betas(3+k+1) = I;
			betas2(3+k+1) = I2;
			s = [s '\t%g'];
			k = k+1;
		end	
		s = [s '\n'];
		fprintf(fid, s, betas);
		fprintf(fid, s, betas2);
	end
	fclose(fid);	
end	