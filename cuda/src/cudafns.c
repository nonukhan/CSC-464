


// number of threads must be the power of 2 closest to but not exceeding n
__global__ void reduce(double* x, double* v, int n) {
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


int find_closest_pow_two(int n) {
	int next = 1;
	int old = 1;

	while (next < n) {
		old = next;
		next *= 2;
	}
	
	return old;
}




double evaluate(node *p, boolean saveit)
{
#ifdef CALLCOUNT
	printf("evaluate\n");
#endif
#ifdef TIMINGS
	get_start_time(eval_tk);
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
#ifdef TIMINGS
		get_stop_time(eval_tk);
#endif
		return sum;
	}
	if(which <= shimotrees) {
		l0gl[which - 1] = sum;
	}
	if (which == 1) {
		maxwhich = 1;
		maxlogl = sum;
#ifdef TIMINGS
		get_stop_time(eval_tk);
#endif
		return sum;
	}
	if (sum > maxlogl) {
		maxwhich = which;
		maxlogl = sum;
	}
#ifdef TIMINGS
	get_stop_time(eval_tk);
#endif
	
	return sum;
}  /* evaluate */




void slopecurv(node *p,double y,double *like,double *slope,double *curve)
{
#ifdef CALLCOUNT
	printf("slopecurve\n");
#endif
#ifdef TIMINGS
	get_start_time(slopecurv_tk);
#endif
	
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

	// rcategs and categs are 1 by default therefore nothing to make parallel :(
	for (i = 0; i < rcategs; i++) {
		for (j = 0; j < categs; j++) {
			tbl[i][j]->orig_zz = exp(tbl[i][j]->rat * lz);
			tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
		}
	}

	// bigger than 1 but less than the # of basepairs
	for (i = 0; i < endsite; i++) {
		k = category[alias[i]-1] - 1;
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
		
		// rcategs = 1
		sumterm = 0.0;
		for (j = 0; j < rcategs; j++) {
			sumterm += probcat[j] * term[i][j];
		}
		lterm = log(sumterm) + p->underflows[i] + q->underflows[i];
		for (j = 0; j < rcategs; j++) {
			term[i][j] = term[i][j] / sumterm;
			slopeterm[i][j] = slopeterm[i][j] / sumterm;
			curveterm[i][j] = curveterm[i][j] / sumterm; 
		}
		sum += aliasweight[i] * lterm;
	}
	// rcategs = 1
	for (i = 0; i < rcategs; i++) {
		thelike[i] = 1.0;
		theslope[i] = 0.0;
		thecurve[i] = 0.0;
	}

	// sites = # of basepairs in seq
	for (i = 0; i < sites; i++) {
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
	
#ifdef TIMINGS
	get_stop_time(slopecurv_tk);
#endif
} /* slopecurv */

