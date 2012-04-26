
double evaluate(node *p, boolean saveit)
{
#ifdef CALLCOUNT
	printf("evaluate\n");
#endif
	
	contribarr tterm;
	double sum, sum2, sumc, y, lz, y1, z1zz, z1yy, prod12, prod1, prod2, prod3, sumterm, lterm;
	long i, j, k, lai;
	node *q;
	sitelike x1, x2;
	
	sum = 0.0;
	q = p->back;
	if ( p->initialized  == false && p->tip == false)  nuview(p);
	if ( q->initialized  == false && q->tip == false)  nuview(q);
	y = p->v;
	lz = -y;

	// categs = rcategs = 1
	for (i = 0; i < rcategs; i++) {
		for (j = 0; j < categs; j++) {
			tbl[i][j]->orig_zz = exp(tbl[i][j]->ratxi * lz);
			tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
			tbl[i][j]->z1zz = tbl[i][j]->z1 * tbl[i][j]->orig_zz;
			tbl[i][j]->z1yy = tbl[i][j]->z1 - tbl[i][j]->z1zz;
		}
	}

	// endsites > 1 but less than sites
	for (i = 0; i < endsite; i++) {
		k = category[alias[i]-1] - 1;

		for (j = 0; j < rcategs; j++) {
			if (y > 0.0) {
				y1 = 1.0 - tbl[j][k]->z1;
				z1zz = tbl[j][k]->z1zz;
				z1yy = tbl[j][k]->z1yy;
			} else {
				y1 = 0.0;
				z1zz = 1.0;
				z1yy = 0.0;
			}
			memcpy(x1, p->x[i][j], sizeof(sitelike));
			prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
				freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
			memcpy(x2, q->x[i][j], sizeof(sitelike));
			prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
				freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
			prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
				(x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
				(x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
				(x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
			prod12 = freqa * x1[0] * x2[0] +
				freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
				freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
				freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
			tterm[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2;
		}

		// rcategs = 1
		sumterm = 0.0;
		for (j = 0; j < rcategs; j++) {
			sumterm += probcat[j] * tterm[j];
		}

		lterm = log(sumterm) + p->underflows[i] + q->underflows[i];

		// rcategs = 1
		for (j = 0; j < rcategs; j++) {
			clai[j] = tterm[j] / sumterm;
		}

		memcpy(contribution[i], clai, rcategs*sizeof(double));
		if (saveit && !auto_ && usertree && (which <= shimotrees)) {
			l0gf[which - 1][i] = lterm;
		}
		sum += aliasweight[i] * lterm;
	}

	// rcategs = 1
	for (j = 0; j < rcategs; j++) {
		like[j] = 1.0;
	}

	// sites = # bp's per seq
	for (i = 0; i < sites; i++) {
		sumc = 0.0;
		for (k = 0; k < rcategs; k++) {
			sumc += probcat[k] * like[k];
		}
		sumc *= lambda;

		if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
			lai = location[ally[i] - 1];
			memcpy(clai, contribution[lai - 1], rcategs*sizeof(double));
			for (j = 0; j < rcategs; j++) {
				nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
			}
		} else {
			for (j = 0; j < rcategs; j++) {
				nulike[j] = ((1.0 - lambda) * like[j] + sumc);
			}
		}
		memcpy(like, nulike, rcategs*sizeof(double));
	}


	// rcategs = 1
	sum2 = 0.0;
	for (i = 0; i < rcategs; i++) {
		sum2 += probcat[i] * like[i];
	}

	sum += log(sum2);
	curtree.likelihood = sum;
	if (!saveit || auto_ || !usertree) {
		return sum;
	}
	if(which <= shimotrees) {
		l0gl[which - 1] = sum;
	}
	if (which == 1) {
		maxwhich = 1;
		maxlogl = sum;
		return sum;
	}
	if (sum > maxlogl) {
		maxwhich = which;
		maxlogl = sum;
	}
	
	return sum;
}  /* evaluate */

/******************************CUDA FN'S********************************/

/*
typedef struct valrec {
	double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1, xiy1xv;
	double *ww, *zz, *wwzz, *vvzz; 
} valrec;
 */

struct _gpu_tbl {
	double *gpu_rat;
	double *gpu_ratxi;
	double *gpu_ratxv;
	double *gpu_orig_zz;
	double *gpu_z1;
	double *gpu_y1;
	double *gpu_z1zz;
	double *gpu_z1yy;
	double *gpu_xiz1;
	double *gpu_xiy1xv;
	double **gpu_ww;
	double **gpu_zz;
	double **gpu_wwzz;
	double **gpu_vvzz;
};

struct _gpu_tbl *gpu_tbl;

// will break if categ != 1 !!
void init_gpu_tbl() {

	gpu_tbl = (struct _gpu_tbl *) malloc(sizeof(struct _gpu_tbl));
	
	cudaMalloc( (void **) &gpu_tbl->gpu_rat, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_ratxi, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_ratxv, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_orig_zz, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_z1, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_y1, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_z1zz, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_z1yy, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_xiz1, rcategs*sizeof(double) );
	cudaMalloc( (void **) &gpu_tbl->gpu_xiy1xv, rcategs*sizeof(double) );
	
}


// tbl_origzz_z1
//
__global__ void _tbl_origzz_z1(double* orig_zz, double* z1,
							   double* rat, double* ratxv, double lz) {
	int tid;
	tid = threadIdx.x;
	
	orig_zz[tid] = exp(rat[tid] * lz);
	z1[tid] = exp(ratxv[tid] * lz);
}

void tbl_origzz_z1(double lz) {
	double *tmp;
	int i, arr_size;
	
	arr_size = rcategs*sizeof(double);
	tmp = (double *) malloc(arr_size);
	
	// copy to device
	for(i=0; i<rcategs; i++) {
		tmp[i] = tbl[i][0]->orig_zz;
	}
	cudaMemcpy(gpu_tbl->gpu_orig_zz, tmp, arr_size, cudaMemcpyHostToDevice);
	for(i=0; i<rcategs; i++) {
		tmp[i] = tbl[i][0]->z1;
	}
	cudaMemcpy(gpu_tbl->gpu_z1, tmp, arr_size, cudaMemcpyHostToDevice);
	for(i=0; i<rcategs; i++) {
		tmp[i] = tbl[i][0]->rat;
	}
	cudaMemcpy(gpu_tbl->gpu_rat, tmp, arr_size, cudaMemcpyHostToDevice);
	for(i=0; i<rcategs; i++) {
		tmp[i] = tbl[i][0]->ratxv;
	}
	cudaMemcpy(gpu_tbl->gpu_ratxv, tmp, arr_size, cudaMemcpyHostToDevice);
	
	// call kernel
	_tbl_origzz_z1<<<1, rcategs>>>(gpu_tbl->gpu_orig_zz, gpu_tbl->gpu_z1,
								   gpu_tbl->gpu_rat, gpu_tbl->gpu_ratxv, lz);
	
	// copy stuff back to host
	cudaMemcpy(tmp, gpu_tbl->gpu_orig_zz, arr_size, cudaMemcpyDeviceToHost);
	for(i=0; i<rcategs; i++) {
		tbl[i][0]->orig_zz = tmp[i];
	}
	cudaMemcpy(tmp, gpu_tbl->gpu_z1, arr_size, cudaMemcpyDeviceToHost);
	for(i=0; i<rcategs; i++) {
		tbl[i][0]->z1 = tmp[i];
	}
	cudaMemcpy(tmp, gpu_tbl->gpu_rat, arr_size, cudaMemcpyDeviceToHost);
	for(i=0; i<rcategs; i++) {
		tbl[i][0]->rat = tmp[i];
	}
	cudaMemcpy(tmp, gpu_tbl->gpu_ratxv, arr_size, cudaMemcpyDeviceToHost);
	for(i=0; i<rcategs; i++) {
		tbl[i][0]->ratxv = tmp[i];
	}
	
	free(tmp);
}
// tbl_origzz_z1


// probcat_sumterm
//
int cuda_pow_two = 0;

int find_closest_pow_two(int n) {
	int next = 1;
	int old = 1;
	
	while (next < n) {
		old = next;
		next *= 2;
	}
	
	return old;
}

// number of threads must be the power of 2 closest to but not exceeding n
__global__ void _prod_reduce(double* x, double* v, int n) {
	int tid, i, ind, num_threads;
	
	num_threads = blockDim.x;
	tid = threadIdx.x;
	
	v[tid] = x[tid] * v[tid];
	// n - num_threads will be leftovers of power of 2
	// add them so we can consider the rest of the reduction as a power of 2
	if(tid < n - num_threads) {
		v[tid] += x[tid + num_threads] * v[tid + num_threads];
	}
	
	__syncthreads();
	
	for(i=1; i < blockDim.x; i *= 2) {
		ind = 2 * i * tid;
		
		if (ind < blockDim.x) {
		    v[ind] += v[ind + i];
		}
		__syncthreads();
	}
}

// used by sumterm and calcterm
double *gpu_probcat;

double *gpu_term;

double probcat_sumterm(double* probcat, double* term) {
	double sum;
	int arr_size;
	
	if (!cuda_pow_two) {
		cuda_pow_two = find_closest_pow_two(rcategs);
	}
	
	arr_size = rcategs*sizeof(double);
	
	if (!gpu_term) {
		cudaMalloc((void **) &gpu_term, arr_size);
	}
	cudaMemcpy(gpu_term, term, arr_size, cudaMemcpyHostToDevice);
	
	if (!gpu_probcat) {
		cudaMalloc((void **) &gpu_probcat, arr_size);
		cudaMemcpy(gpu_probcat, probcat, arr_size, cudaMemcpyHostToDevice);
	}
	
	_prod_reduce<<<1, cuda_pow_two>>>(gpu_term, gpu_probcat, rcategs);
	
	cudaMemcpy(&sum, gpu_term, sizeof(double), cudaMemcpyDeviceToHost);
	
	return sum;
}
// probcat_sumterm



// calc_terms
//
__global__ void _calc_terms(double *term, double *slopeterm, 
							double *curveterm, double sumterm) {
	int i;
	i = threadIdx.x;
	
	term[i] = term[i] / sumterm;
	slopeterm[i] = slopeterm[i] / sumterm;
	curveterm[i] = curveterm[i] / sumterm; 
}


double *gpu_slopeterm, *gpu_curveterm;

void calc_terms(double *termm, double *slopetermm, 
				double *curvetermm, double sumterm) {
	int arr_size;
	
	arr_size = rcategs*sizeof(double);
	
	if (!gpu_term) {
		cudaMalloc((void **) &gpu_term, arr_size);
	}
	cudaMemcpy(gpu_term, termm, arr_size, cudaMemcpyHostToDevice);
	
	if (!gpu_slopeterm) {
		cudaMalloc((void **) &gpu_slopeterm, arr_size);
	}
	cudaMemcpy(gpu_slopeterm, slopetermm, arr_size, cudaMemcpyHostToDevice);
	
	if (!gpu_curveterm) {
		cudaMalloc((void **) &gpu_curveterm, arr_size);
	}
	cudaMemcpy(gpu_curveterm, curvetermm, arr_size, cudaMemcpyHostToDevice);

	_calc_terms<<<1, rcategs>>>(gpu_term, gpu_slopeterm, gpu_curveterm, sumterm);
	
	cudaMemcpy(termm, gpu_term, arr_size, cudaMemcpyDeviceToHost);
	cudaMemcpy(slopetermm, gpu_slopeterm, arr_size, cudaMemcpyDeviceToHost);
	cudaMemcpy(curvetermm, gpu_curveterm, arr_size, cudaMemcpyDeviceToHost);
}
// calc_terms


// calc_sums
//

// number of threads must be the power of 2 closest to but not exceeding n
__global__ void _calc_sums(double* like, double* slope, double* curve, double* probcat, int n) {
	int tid, i, ind, num_threads;
	
	num_threads = blockDim.x;
	tid = threadIdx.x;
	
	like[tid] = like[tid] * probcat[tid];
	slope[tid] = slope[tid] * probcat[tid];
	curve[tid] = curve[tid] * probcat[tid];

	// n - num_threads will be leftovers of power of 2
	// add them so we can consider the rest of the reduction as a power of 2
	if(tid < n - num_threads) {		
		like[tid] += like[tid + num_threads] * probcat[tid + num_threads];
		slope[tid] += slope[tid + num_threads] * probcat[tid + num_threads];
		curve[tid] += curve[tid + num_threads] * probcat[tid + num_threads];
	}
	
	__syncthreads();
	
	for(i=1; i < blockDim.x; i *= 2) {
		ind = 2 * i * tid;
		
		if (ind < blockDim.x) {			
			like[tid] += like[tid + i];
			slope[tid] += slope[tid + i];
			curve[tid] += curve[tid + i];
		}
		__syncthreads();
	}
}

double *gpu_like, *gpu_slope, *gpu_curve;
void calc_sums(double *sumc, double *sumcs, double *sumcc, double *probcat, 
			   double *thelike, double *theslope, double *thecurve, double lambda) {
	int arr_size;
	double *gpu_like, *gpu_slope, *gpu_curve;
	
	arr_size = rcategs*sizeof(double);
	
	if (!gpu_probcat) {
		cudaMalloc((void **) &gpu_probcat, arr_size);
		cudaMemcpy(gpu_probcat, probcat, arr_size, cudaMemcpyHostToDevice);
	}
	
	if (!gpu_like) {
		cudaMalloc((void **) &gpu_like, arr_size);
	}
	cudaMemcpy(gpu_like, thelike, arr_size, cudaMemcpyHostToDevice);
	
	if (!gpu_slope) {
		cudaMalloc((void **) &gpu_slope, arr_size);
	}
	cudaMemcpy(gpu_slope, theslope, arr_size, cudaMemcpyHostToDevice);
	
	if (!gpu_curve) {
		cudaMalloc((void **) &gpu_curve, arr_size);
	}
	cudaMemcpy(gpu_curve, thecurve, arr_size, cudaMemcpyHostToDevice);
	
	_calc_sums<<<1, cuda_pow_two>>>(gpu_slope, gpu_like, gpu_curve, gpu_probcat, rcategs);
	
	cudaMemcpy(sumc, gpu_like, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(sumcs, gpu_slope, sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(sumcc, gpu_curve, sizeof(double), cudaMemcpyDeviceToHost);

	*sumc *= lambda;
	*sumcs *= lambda;
	*sumcc *= lambda;
}

void gpu_free() {
	
	cudaFree(gpu_like);
	cudaFree(gpu_slope);
	cudaFree(gpu_curve);
	cudaFree(gpu_term);
	cudaFree(gpu_slopeterm);
	cudaFree(gpu_curveterm);
	cudaFree(gpu_probcat);
}

//

/******************************CUDA FN'S********************************/

void slopecurv(node *p,double y,double *like,double *slope,double *curve)
{
#ifdef CALLCOUNT
	printf("slopecurve\n");
#endif
#ifdef SCTIMINGS
	get_start_time(slopecurv_tk);
#endif
	
	if (!gpu_tbl) {
		init_gpu_tbl();
	}
	//update_gpu_tbl();
	
	/* compute log likelihood, slope and curvature at node p */
	long i, j, k, lai;
	double sum, sumc, sumterm, lterm, sumcs, sumcc, sum2, slope2, curve2,
	temp;
	double lz, zz, z1, zzs, z1s, zzc, z1c, aa, bb, cc,
	prod1, prod2, prod12, prod3;
	contribarr thelike, nulike, nuslope, nucurve,
    theslope, thecurve, clai, cslai, cclai;
	node *q;
	sitelike x1, x2;
	
	
	q = p->back;
	sum = 0.0;
	lz = -y;

	
	/*
	for (i = 0; i < rcategs; i++) {
		for (j = 0; j < categs; j++) {
			tbl[i][j]->orig_zz = exp(tbl[i][j]->rat * lz);
			tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
		}
	}
	 */
	//printf("tbl_origzz_z1\n");
	tbl_origzz_z1(lz);
	
	// bigger than 1 but less than the # of basepairs
	for (i = 0; i < endsite; i++) {
		k = category[alias[i]-1] - 1;
		
		//
		for (j = 0; j < rcategs; j++) {
			if (y > 0.0) {
				zz = tbl[j][k]->orig_zz;
				z1 = tbl[j][k]->z1;
			} else {
				zz = 1.0;
				z1 = 1.0;
			}
			zzs = -tbl[j][k]->rat * zz ;
			z1s = -tbl[j][k]->ratxv * z1 ;
			temp = tbl[j][k]->rat;
			zzc = temp * temp * zz;
			temp = tbl[j][k]->ratxv;
			z1c = temp * temp * z1;
			memcpy(x1, p->x[i][j], sizeof(sitelike));
			prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
           		freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
			memcpy(x2, q->x[i][j], sizeof(sitelike));
			prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
           		freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
			prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
				(x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
				(x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
				(x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
			prod12 = freqa * x1[0] * x2[0] +
				freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
				freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
				freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
			aa = prod12 - prod3;
			bb = prod3 - prod1*prod2;
			cc = prod1 * prod2;
			term[i][j] = zz * aa + z1 * bb + cc;
			slopeterm[i][j] = zzs * aa + z1s * bb;
			curveterm[i][j] = zzc * aa + z1c * bb;
		}
		//
		
		/*
		sumterm = 0.0;
		for (j = 0; j < rcategs; j++) {
			sumterm += probcat[j] * term[i][j];
		}
		 */
		//printf("probcat_sumterm\n");
		sumterm = probcat_sumterm(probcat, term[i]);

		lterm = log(sumterm) + p->underflows[i] + q->underflows[i];
		
		/*
		for (j = 0; j < rcategs; j++) {
			term[i][j] = term[i][j] / sumterm;
			slopeterm[i][j] = slopeterm[i][j] / sumterm;
			curveterm[i][j] = curveterm[i][j] / sumterm; 
		}
		 */	
		//printf("calc_terms\n");
		calc_terms(term[i], slopeterm[i], curveterm[i], sumterm);
		
		sum += aliasweight[i] * lterm;
	}
	
	// one memcpy is already slower
	for (i = 0; i < rcategs; i++) {
		thelike[i] = 1.0;
		theslope[i] = 0.0;
		thecurve[i] = 0.0;
	}

	// sites = # of basepairs in seq
	for (i = 0; i < sites; i++) {
		/*
		sumc = 0.0;
		sumcs = 0.0;
		sumcc = 0.0;
		for (k = 0; k < rcategs; k++) {
			sumc += probcat[k] * thelike[k];
			sumcs += probcat[k] * theslope[k];
			sumcc += probcat[k] * thecurve[k];
		}
		sumc *= lambda;
		sumcs *= lambda;
		sumcc *= lambda;
		 */
		calc_sums(&sumc, &sumcs, &sumcc, probcat, thelike, theslope, thecurve, lambda);
		
		if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
			lai = location[ally[i] - 1];
			memcpy(clai, term[lai - 1], rcategs*sizeof(double));
			memcpy(cslai, slopeterm[lai - 1], rcategs*sizeof(double));
			memcpy(cclai, curveterm[lai - 1], rcategs*sizeof(double));
			if (weight[i] > 1) {
				
				for (j = 0; j < rcategs; j++) {
					if (clai[j] > 0.0) {
						clai[j] = exp(weight[i]*log(clai[j]));
					} else {
						clai[j] = 0.0;
					}
					if (cslai[j] > 0.0) {
						cslai[j] = exp(weight[i]*log(cslai[j]));
					} else {
						cslai[j] = 0.0;
					}
					if (cclai[j] > 0.0) {
						cclai[j] = exp(weight[i]*log(cclai[j])); 
					} else {
						cclai[j] = 0.0;
					}
				}
				
			}
			for (j = 0; j < rcategs; j++) {
				nulike[j] = ((1.0 - lambda) * thelike[j] + sumc) * clai[j];
				nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs) * clai[j]
					+ ((1.0 - lambda) * thelike[j] + sumc) * cslai[j];
				nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc) * clai[j]
					+ 2.0 * ((1.0 - lambda) * theslope[j] + sumcs) * cslai[j]
					+ ((1.0 - lambda) * thelike[j] + sumc) * cclai[j];
			}
		} else {
			for (j = 0; j < rcategs; j++) {
				nulike[j] = ((1.0 - lambda) * thelike[j] + sumc);
				nuslope[j] = ((1.0 - lambda) * theslope[j] + sumcs);
				nucurve[j] = ((1.0 - lambda) * thecurve[j] + sumcc);
			}
		}
		// changes thelike, etc.
		memcpy(thelike, nulike, rcategs*sizeof(double));
		memcpy(theslope, nuslope, rcategs*sizeof(double));
		memcpy(thecurve, nucurve, rcategs*sizeof(double));
	}
	sum2 = 0.0;
	slope2 = 0.0;
	curve2 = 0.0;

	// rcategs is 1, no parallel :(
	for (i = 0; i < rcategs; i++) {
		sum2 += probcat[i] * thelike[i];
		slope2 += probcat[i] * theslope[i];
		curve2 += probcat[i] * thecurve[i];
	}

	sum += log(sum2);
	(*like) = sum;
	(*slope) = slope2 / sum2;
	
	/* Expressed in terms of *slope to prevent overflow */
	(*curve) = curve2 / sum2 - *slope * *slope;
	
#ifdef SCTIMINGS
	get_stop_time(slopecurv_tk);
#endif
} /* slopecurv */

