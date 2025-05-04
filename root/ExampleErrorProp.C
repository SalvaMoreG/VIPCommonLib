
int Main()
{
	// step 1
	// 
	float mean_100pc_set_1 = 45.2;
	float sigma_100pc_set_1 = 0.195;
	float mean_100pc_set_2 = 10;
	float sigma_100pc_set_2 = 0.01;

	// float f = 4.52; 
	// float sigma_f = 0.02; 
	float f = mean_100pc_set_1/mean_100pc_set_2;

	float tmp1 = sigma_100pc_set_1/mean_100pc_set_1;
	float tmp2 = sigma_100pc_set_2/mean_100pc_set_2;
	float tmp1_sq = tmp1 * tmp1;
	float tmp2_sq = tmp2 * tmp2;

	float sigma_f = f * sqrt(tmp1_sq + tmp2_sq);
	cout << "mean 100 % set 1: " << mean_100pc_set_1 << " +- " << sigma_100pc_set_1 << endl;
	cout << "mean 100 % set 2: " << mean_100pc_set_2 << " +- " << sigma_100pc_set_2 << endl;
	cout << "normalization factor: " << f << " +- " << sigma_f << endl;
	cout << endl;

	// step 2
	// 
	float mean_old = 2.0; 
	float sigma_old = 0.2;

	float mean_normalized = mean_old * f;
	tmp1 = sigma_f/f;
	tmp2 = sigma_old/mean_old;
	tmp1_sq = tmp1 * tmp1;
	tmp2_sq = tmp2 * tmp2;

	float sigma_normalized = mean_normalized * sqrt(tmp1_sq + tmp2_sq);
	cout << "normalization factor: " << f << " +- " << sigma_f << endl;
	cout << "old value: " << mean_old << " +- " << sigma_old << endl;
	cout << "normalized value: " << mean_normalized << " +- " << sigma_normalized << endl;

	return 0;

}

