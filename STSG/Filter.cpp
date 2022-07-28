#include "Filter.h"

#include <math.h>

#include <algorithm>
using namespace std;

void STSG_filter(int win, float sampcorr, float *img_NDVI, float *img_QA, float *reference_data, int ns, int nl, int nb, int ny, int snow_address, float *vector_out) {
    for (int i = win; i < ns - win; i++) {
        for (int j = win; j < nl - win; j++) {
            int samp = 0;
            int aap = 0;
            int *similar_index = new int[2 * (2 * win + 1) * (2 * win + 1)];
            float *slope_intercept = new float[2 * (2 * win + 1) * (2 * win + 1)];
            float *corr_similar = new float[(2 * win + 1) * (2 * win + 1)];
            for (int y = 0; y < ny; y++) {
                float *VI_raw = new float[nb];
                float *vector_QA = new float[nb];
                float *res_VI_raw = new float[nb];
                for (int k = 0; k < nb; k++) {
                    int ind = gdis = i + j * ns + k * ns * nl + y * ns * nl * nb;
                    VI_raw[k] = img_NDVI[ind];
                    vector_QA[k] = img_QA[ind];
                    res_VI_raw[k] = img_NDVI[ind];
                }

                for (int m = 0; m < nb - 1; m++) {
                    for (int n = m + 1; n < nb; n++) {
                        if (res_VI_raw[m] < res_VI_raw[n]) {
                            float temp = res_VI_raw[m];
                            res_VI_raw[m] = res_VI_raw[n];
                            res_VI_raw[n] = temp;
                        }
                    }
                }
                if (((res_VI_raw[0] + res_VI_raw[1] + res_VI_raw[2]) / 3) > 0.15)  //若小于等于呢？
                {
                    int indic = 0;
                    // searching similar pixels
                    if (y == 0) {
                        float *corr_res = new float[(2 * win + 1) * (2 * win + 1)];
                        float *Intercept_res = new float[(2 * win + 1) * (2 * win + 1)];
                        float *Slope_res = new float[(2 * win + 1) * (2 * win + 1)];
                        float *new_corr_similar_res = new float[(2 * win + 1) * (2 * win + 1)];
                        int count_corr_Slope = 0;
                        for (int si = -1 * win; si <= win; si++) {
                            for (int sj = -1 * win; sj <= win; sj++) {
                                if (reference_data[i + si + (j + sj) * ns + 3 * ns * nl] != 0)  //为什么是3？不是怎么办？
                                {
                                    double x_sum = 0;
                                    double y_sum = 0;
                                    double x_mean = 0;
                                    double y_mean = 0;
                                    double xy_sum = 0;
                                    double x2_sum = 0;
                                    double y2_sum = 0;
                                    for (int k = 0; k < nb; k++) {
                                        x_sum += reference_data[i + si + (j + sj) * ns + k * ns * nl];
                                        y_sum += reference_data[i + j * ns + k * ns * nl];
                                    }
                                    x_mean = x_sum / nb;
                                    y_mean = y_sum / nb;
                                    for (int k = 0; k < nb; k++) {
                                        int ind = i + si + (j + sj) * ns + k * ns * nl;
                                        xy_sum += (reference_data[ind] - x_mean) * (reference_data[i + j * ns + k * ns * nl] - y_mean);
                                        x2_sum += (reference_data[ind] - x_mean) * (reference_data[ind] - x_mean);
                                        y2_sum += (reference_data[i + j * ns + k * ns * nl] - y_mean) * (reference_data[i + j * ns + k * ns * nl] - y_mean);
                                    }
                                    corr_res[si + win + (sj + win) * (2 * win + 1)] = xy_sum / sqrt(x2_sum * y2_sum);

                                    x_sum = 0;
                                    y_sum = 0;
                                    x_mean = 0;
                                    y_mean = 0;
                                    xy_sum = 0;
                                    x2_sum = 0;
                                    y2_sum = 0;
                                    int count_tempQA = 0;
                                    for (int k = 0; k < nb; k++) {
                                        for (int yeari = 0; yeari < ny; yeari++) {
                                            int ind = i + si + (j + sj) * ns + k * ns * nl; 
                                            if (img_QA[i + j * ns + k * ns * nl + yeari * ns * nl * nb] <= 1 && img_QA[ind + yeari * ns * nl * nb] <= 1) {
                                                if (reference_data[i + j * ns + k * ns * nl] != 0 && reference_data[ind] != 0) {
                                                    x_sum += img_NDVI[ind + yeari * ns * nl * nb] / reference_data[ind];
                                                    y_sum += img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb] / reference_data[i + j * ns + k * ns * nl];
                                                    xy_sum += (img_NDVI[ind + yeari * ns * nl * nb] / reference_data[ind]) * (img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb] / reference_data[i + j * ns + k * ns * nl]);
                                                    x2_sum += (img_NDVI[ind + yeari * ns * nl * nb] / reference_data[ind]) * (img_NDVI[ind + yeari * ns * nl * nb] / reference_data[ind]);
                                                    y2_sum += (img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb] / reference_data[i + j * ns + k * ns * nl]) * (img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb] / reference_data[i + j * ns + k * ns * nl]);
                                                    count_tempQA++;
                                                }
                                            }
                                        }
                                    }
                                    if (count_tempQA >= 30) {
                                        Intercept_res[si + win + (sj + win) * (2 * win + 1)] = (x2_sum * y_sum - x_sum * xy_sum) / (x2_sum * count_tempQA - x_sum * x_sum);
                                        Slope_res[si + win + (sj + win) * (2 * win + 1)] = (xy_sum * count_tempQA - x_sum * y_sum) / (x2_sum * count_tempQA - x_sum * x_sum);

                                        x_mean = x_sum / count_tempQA;
                                        y_mean = y_sum / count_tempQA;
                                        xy_sum = 0;
                                        x2_sum = 0;
                                        y2_sum = 0;
                                        for (int k = 0; k < nb; k++) {
                                            for (int yeari = 0; yeari < ny; yeari++) {
                                                int ind = i + si + (j + sj) * ns + k * ns * nl;

                                                if (img_QA[i + j * ns + k * ns * nl + yeari * ns * nl * nb] <= 1 && img_QA[ind + yeari * ns * nl * nb] <= 1) {

                                                    if (reference_data[i + j * ns + k * ns * nl] != 0 && reference_data[ind] != 0) {
                                                        xy_sum += (img_NDVI[ind + yeari * ns * nl * nb] / reference_data[ind] - x_mean) * (img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb] / reference_data[i + j * ns + k * ns * nl] - y_mean);
                                                        x2_sum += (img_NDVI[ind + yeari * ns * nl * nb] / reference_data[ind] - x_mean) * (img_NDVI[ind + yeari * ns * nl * nb] / reference_data[ind] - x_mean);
                                                        y2_sum += (img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb] / reference_data[i + j * ns + k * ns * nl] - y_mean) * (img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb] / reference_data[i + j * ns + k * ns * nl] - y_mean);
                                                    }
                                                }
                                            }
                                        }
                                        new_corr_similar_res[si + win + (sj + win) * (2 * win + 1)] = xy_sum / sqrt(x2_sum * y2_sum);
                                    } else {
                                        Slope_res[si + win + (sj + win) * (2 * win + 1)] = 0;
                                        Intercept_res[si + win + (sj + win) * (2 * win + 1)] = 0;
                                        new_corr_similar_res[si + win + (sj + win) * (2 * win + 1)] = 0;
                                    }
                                } else {
                                    corr_res[si + win + (sj + win) * (2 * win + 1)] = 0;
                                    Slope_res[si + win + (sj + win) * (2 * win + 1)] = 0;
                                    Intercept_res[si + win + (sj + win) * (2 * win + 1)] = 0;
                                    new_corr_similar_res[si + win + (sj + win) * (2 * win + 1)] = 0;
                                }
                            }
                        }

                        for (int m = 0; m < (2 * win + 1) * (2 * win + 1); m++) {
                            if (corr_res[m] >= sampcorr && Slope_res[m] != 0.)
                                count_corr_Slope++;
                        }
                        if (count_corr_Slope >= 2) {
                            int *new_corr = new int[(2 * win + 1) * (2 * win + 1)];
                            for (int m = 0; m < (2 * win + 1) * (2 * win + 1); m++)
                                new_corr[m] = m;
                            for (int m = 0; m < (2 * win + 1) * (2 * win + 1) - 1; m++) {
                                for (int n = m + 1; n < (2 * win + 1) * (2 * win + 1); n++) {
                                    if (corr_res[new_corr[m]] < corr_res[new_corr[n]]) {
                                        int temp = new_corr[m];
                                        new_corr[m] = new_corr[n];
                                        new_corr[n] = temp;
                                    }
                                }
                            }

                            samp = count_corr_Slope - 1;
                            for (int m = 0; m < samp; m++) {
                                similar_index[1 + m * 2] = int(new_corr[m + 1] / (2 * win + 1)) + j - win;
                                similar_index[m * 2] = new_corr[m + 1] - int(new_corr[m + 1] / (2 * win + 1)) * (2 * win + 1) + i - win;
                                slope_intercept[1 + m * 2] = Slope_res[int(new_corr[m + 1] / (2 * win + 1)) + (2 * win + 1) * (new_corr[m + 1] - int(new_corr[m + 1] / (2 * win + 1)) * (2 * win + 1))];
                                slope_intercept[m * 2] = Intercept_res[int(new_corr[m + 1] / (2 * win + 1)) + (2 * win + 1) * (new_corr[m + 1] - int(new_corr[m + 1] / (2 * win + 1)) * (2 * win + 1))];
                                corr_similar[m] = new_corr_similar_res[int(new_corr[m + 1] / (2 * win + 1)) + (2 * win + 1) * (new_corr[m + 1] - int(new_corr[m + 1] / (2 * win + 1)) * (2 * win + 1))];
                            }
                            aap = 1;

                            delete[] new_corr;
                        } else
                            aap = 0;

                        delete[] corr_res;
                        delete[] Intercept_res;
                        delete[] Slope_res;
                        delete[] new_corr_similar_res;
                    }

                    // generate the trend curve
                    float *VI_init = new float[nb];
                    if (aap == 1) {
                        for (int k = 0; k < nb; k++) {
                            float *temp_NDVI = new float[samp];
                            int count_temp_NDVI = 0;
                            float total_new_corr_similar = 0;
                            float total_new_temp = 0;
                            for (int m = 0; m < samp; m++) {
                                if (img_QA[similar_index[m * 2] + similar_index[1 + m * 2] * ns + k * ns * nl + y * ns * nl * nb] <= 1 && reference_data[similar_index[m * 2] + similar_index[1 + m * 2] * ns + k * ns * nl] != 0) {
                                    float new_ratio = img_NDVI[similar_index[m * 2] + similar_index[1 + m * 2] * ns + k * ns * nl + y * ns * nl * nb] / reference_data[similar_index[m * 2] + similar_index[1 + m * 2] * ns + k * ns * nl];
                                    temp_NDVI[m] = (slope_intercept[m * 2] + new_ratio * slope_intercept[1 + m * 2]) * reference_data[i + j * ns + k * ns * nl];
                                    if (temp_NDVI[m] >= 1 || temp_NDVI[m] <= -0.2)
                                        temp_NDVI[m] = 0;
                                } else
                                    temp_NDVI[m] = 0;

                                if (temp_NDVI[m] != 0) {
                                    count_temp_NDVI++;
                                    total_new_corr_similar += corr_similar[m];
                                }
                            }
                            if (count_temp_NDVI != 0) {
                                for (int m = 0; m < samp; m++) {
                                    if (temp_NDVI[m] != 0)
                                        total_new_temp += corr_similar[m] / total_new_corr_similar * temp_NDVI[m];
                                }
                                VI_init[k] = total_new_temp;
                            } else
                                VI_init[k] = 0;

                            delete[] temp_NDVI;
                        }

                        // generating the VI_init
                        int count_VI_init = 0;
                        int count_VI_init_no = 0;
                        int *res_VI_init = new int[nb];
                        int *res_VI_init_no = new int[nb];
                        int count_conres = 0;
                        for (int k = 0; k < nb; k++) {
                            if (VI_init[k] != 0) {
                                res_VI_init[count_VI_init++] = k;
                                if (count_VI_init > 1 && k - res_VI_init[count_VI_init - 2] >= 3)
                                    count_conres++;
                            } else
                                res_VI_init_no[count_VI_init_no++] = k;
                        }
                        if (count_VI_init >= nb / 2 && count_conres == 0) {
                            for (int m = 0; m < count_VI_init_no; m++) {
                                int start = 0;
                                if (res_VI_init_no[m] < res_VI_init[0])
                                    start = 0;
                                else if (res_VI_init_no[m] > res_VI_init[count_VI_init - 1])
                                    start = count_VI_init - 4;
                                else {
                                    for (int n = 0; n < count_VI_init - 1; n++) {
                                        if (res_VI_init[n] < res_VI_init_no[m] && res_VI_init_no[m] < res_VI_init[n + 1]) {
                                            if (n - 1 < 0)
                                                start = 0;
                                            else if (count_VI_init - n < 4)
                                                start = count_VI_init - 4;
                                            else
                                                start = n - 1;
                                            break;
                                        }
                                    }
                                }

                                float x[4], y[4];
                                for (int n = 0; n < 4; n++) {
                                    x[n] = res_VI_init[start + n];
                                    y[n] = VI_init[res_VI_init[start + n]];
                                }

                                float sig, p;
                                float y2[4] = {0};
                                float u[4] = {0};
                                for (int n = 1; n < 3; n++) {
                                    sig = (x[n] - x[n - 1]) / (x[n + 1] - x[n - 1]);
                                    p = sig * y2[n - 1] + 2;
                                    y2[n] = (sig - 1) / p;
                                    u[n] = (y[n + 1] - y[n]) / (x[n + 1] - x[n]) - (y[n] - y[n - 1]) / (x[n] - x[n - 1]);
                                    u[n] = (6.0 * u[n] / (x[n + 1] - x[n - 1]) - sig * u[n - 1]) / p;
                                }
                                for (int n = 3; n >= 0; n--)
                                    y2[n] = y2[n] * y2[n + 1] + u[n];

                                int klo = 0;
                                int khi = 3;
                                while (khi - klo > 1)  //二分法查找x所在区间段
                                {
                                    int k = (khi + klo) >> 1;
                                    if (x[k] > res_VI_init_no[m])
                                        khi = k;
                                    else
                                        klo = k;
                                }

                                float h = x[khi] - x[klo];
                                float a = (x[khi] - res_VI_init_no[m]) / h;
                                float b = (res_VI_init_no[m] - x[klo]) / h;

                                VI_init[res_VI_init_no[m]] = a * y[klo] + b * y[khi] + ((a * a * a - a) * y2[klo] + (b * b * b - b) * y2[khi]) * h * h / 6.0;
                            }
                            indic = 1;
                        } else
                            indic = 0;

                        delete[] res_VI_init;
                        delete[] res_VI_init_no;
                    } else
                        indic = 0;

                    // STSG
                    if (indic == 1) {
                        // processing contaminated NDVI by snow
                        if (snow_address == 1) {
                            int count_vector_QA = 0;
                            for (int k = 0; k < nb; k++) {
                                if (vector_QA[k] == 2)
                                    count_vector_QA++;
                            }
                            if (count_vector_QA != 0) {
                                int bv_count = 0;
                                float bv_total = 0.;
                                for (int yeari = 0; yeari < ny; yeari++) {
                                    for (int k = 0; k < 6; k++) {
                                        if (img_QA[i + j * ns + k * ns * nl + yeari * ns * nl * nb] <= 1) {
                                            bv_count++;
                                            bv_total += img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb];
                                        }
                                    }
                                }
                                if (bv_count != 0) {
                                    float bv = bv_total / bv_count;
                                    for (int k = 0; k < nb; k++) {
                                        if (vector_QA[k] == 2) {
                                            VI_raw[k] = bv;
                                            VI_init[k] = bv;
                                        }
                                    }
                                }
                            }
                        }

                        // Calculate the weights for each point
                        float gdis = 0.0;
                        float *fl = new float[nb];
                        int count_fl = 0;
                        float mean_fl = 0;
                        for (int k = 0; k < nb; k++) {
                            if (vector_QA[k] == 0 || vector_QA[k] == 1) {
                                fl[k] = VI_raw[k] - VI_init[k];
                                count_fl++;
                                mean_fl += fl[k];
                            } else
                                fl[k] = -1.0;
                        }
                        if (count_fl != 0)
                            mean_fl = mean_fl / count_fl;
                        for (int k = 0; k < nb; k++) {
                            if (fl[k] == -1.0)
                                fl[k] = mean_fl;
                        }

                        for (int k = 0; k < nb; k++) {
                            float min_fl = 0;
                            float max_fl = 0;
                            for (int m = 0; m < nb; m++) {
                                if (max_fl < fl[m])
                                    max_fl = fl[m];
                                if (min_fl > fl[m])
                                    min_fl = fl[m];
                            }
                            fl[k] = (fl[k] - min_fl) / (max_fl - min_fl);
                            if ((VI_raw[k] - VI_init[k]) >= 0)
                                gdis += fl[k] * (VI_raw[k] - VI_init[k]);
                            else
                                gdis += fl[k] * (VI_init[k] - VI_raw[k]);
                        }

                        float *ra4 = new float[nb];
                        float *pre = new float[nb];
                        float ormax = gdis;
                        for (int k = 0; k < nb; k++) {
                            if (vector_QA[k] == 0)
                                VI_init[k] = VI_raw[k];
                            if (vector_QA[k] != 0 && vector_QA[k] != 1)
                                VI_raw[k] = VI_init[k];
                        }

                        int loop_times = 1;
                        while (gdis <= ormax && loop_times < 50) {
                            loop_times += 1;
                            for (int k = 0; k < nb; k++) {
                                ra4[k] = (VI_raw[k] >= VI_init[k]) ? VI_raw[k] : VI_init[k];
                                pre[k] = VI_init[k];
                            }
                            // The Savitzky - Golay fitting
                            // savgolFilter = SAVGOL(4, 4, 0, 6); set the window width(4, 4) and degree(6) for repetition
                            double savgolFilter[] = {-0.00543880, 0.0435097, -0.152289, 0.304585, 0.619267, 0.304585, -0.152289, 0.0435097, -0.00543880};
                            int savgolFilterW = sizeof(savgolFilter) / sizeof(savgolFilter[0]);
                            int ra4W = nb;
                            float *new_ra4 = new float[savgolFilterW + ra4W - 1];
                            for (int m = 0; m < savgolFilterW + ra4W - 1; m++) {
                                if (m < (savgolFilterW - 1) / 2)
                                    new_ra4[m] = ra4[0];
                                else if (m > ((savgolFilterW - 1) / 2 + ra4W - 1))
                                    new_ra4[m] = ra4[ra4W - 1];
                                else
                                    new_ra4[m] = ra4[m - (savgolFilterW - 1) / 2];
                            }
                            for (int m = 0; m < ra4W; m++) {
                                float temp = 0;
                                for (int n = 0; n < savgolFilterW; n++)
                                    temp += savgolFilter[n] * new_ra4[m + n];
                                VI_init[m] = temp;
                            }
                            delete[] new_ra4;
                            ormax = gdis;
                            // Calculate the fitting-effect index
                            gdis = 0.0;
                            for (int k = 0; k < nb; k++) {
                                if ((VI_raw[k] - VI_init[k]) >= 0)
                                    gdis += fl[k] * (VI_raw[k] - VI_init[k]);
                                else
                                    gdis += fl[k] * (VI_init[k] - VI_raw[k]);
                            }
                        }
                        delete[] fl;
                        delete[] ra4;

                        float *vec_fil = new float[nb];
                        for (int k = 0; k < nb; k++)
                            vec_fil[k] = pre[k];
                        for (int smi = 0; smi < nb - 4; smi++) {
                            float a1 = vec_fil[smi];
                            float a2 = vec_fil[smi + 1];
                            float a3 = vec_fil[smi + 2];
                            float a4 = vec_fil[smi + 3];
                            float a5 = vec_fil[smi + 4];
                            if ((a1 > a2) && (a2 < a3) && (a3 > a4) && (a4 < a5)) {
                                pre[smi + 1] = (a1 + a3) / 2.0;
                                pre[smi + 3] = (a3 + a5) / 2.0;
                            }
                        }
                        delete[] vec_fil;

                        for (int k = 0; k < nb; k++)
                            vector_out[i - win + (j - win) * (ns - 2 * win) + k * (ns - 2 * win) * (nl - 2 * win) + y * (ns - 2 * win) * (nl - 2 * win) * nb] = pre[k];
                        delete[] pre;
                    }
                    delete[] VI_init;

                    // SG filter
                    if (indic == 0) {
                        // processing contaminated NDVI by snow
                        if (snow_address == 1) {
                            int count_vector_QA = 0;
                            for (int k = 0; k < nb; k++) {
                                if (vector_QA[k] == 2)
                                    count_vector_QA++;
                            }
                            if (count_vector_QA != 0) {
                                int bv_count = 0;
                                float bv_total = 0.;
                                for (int yeari = 0; yeari < ny; yeari++) {
                                    for (int k = 0; k < 6; k++) {
                                        if (img_QA[i + j * ns + k * ns * nl + yeari * ns * nl * nb] <= 1) {
                                            bv_count++;
                                            bv_total += img_NDVI[i + j * ns + k * ns * nl + yeari * ns * nl * nb];
                                        }
                                    }
                                }
                                if (bv_count != 0) {
                                    float bv = bv_total / bv_count;
                                    for (int k = 0; k < nb; k++) {
                                        if (vector_QA[k] == 2) {
                                            VI_raw[k] = bv;
                                        }
                                    }
                                }
                            }
                        }

                        int count_vector_QA = 0;
                        int *res_vector_QA = new int[nb];
                        for (int k = 0; k < nb; k++) {
                            if (vector_QA[k] <= 2)
                                res_vector_QA[count_vector_QA++] = k;
                        }
                        if (count_vector_QA < nb) {
                            double k;
                            int x = 0;
                            int l = 0;
                            if (res_vector_QA[count_vector_QA - 1] != nb - 1) {
                                k = (VI_raw[res_vector_QA[count_vector_QA - 1]] - VI_raw[res_vector_QA[count_vector_QA - 2]]) / (res_vector_QA[count_vector_QA - 1] - res_vector_QA[count_vector_QA - 2]);
                                VI_raw[nb - 1] = VI_raw[res_vector_QA[count_vector_QA - 1]] + k * (nb - 1 - res_vector_QA[count_vector_QA - 1]);
                                res_vector_QA[count_vector_QA++] = nb - 1;
                            }
                            int count_res_vector_QA = count_vector_QA;
                            while (count_vector_QA < nb && x < count_res_vector_QA - 1) {
                                l = res_vector_QA[x + 1] - res_vector_QA[x];
                                int n = 1;
                                while (l > 1) {
                                    k = (VI_raw[res_vector_QA[x + 1]] - VI_raw[res_vector_QA[x]]) / (res_vector_QA[x + 1] - res_vector_QA[x]);
                                    VI_raw[res_vector_QA[x] + n] = VI_raw[res_vector_QA[x]] + k * n;
                                    count_vector_QA++;
                                    n++;
                                    l--;
                                }
                                x++;
                            }
                            if (res_vector_QA[0] != 0) {
                                k = (VI_raw[res_vector_QA[1]] - VI_raw[res_vector_QA[0]]) / (res_vector_QA[1] - res_vector_QA[0]);
                                l = res_vector_QA[0];
                                int n = 0;
                                do {
                                    VI_raw[n] = VI_raw[res_vector_QA[0]] - k * l;
                                    n++;
                                    l--;
                                } while (l >= 1);
                            }
                        }
                        delete[] res_vector_QA;

                        float *rst = new float[nb];
                        // savgolFilter = SAVGOL(4, 4, 0, 2); set the window width(4, 4) and degree(2) for computing trend curve
                        double savgolFilter[] = {-0.0909091, 0.0606061, 0.168831, 0.233766, 0.255411, 0.233766, 0.168831, 0.0606061, -0.0909091};
                        int savgolFilterW = sizeof(savgolFilter) / sizeof(savgolFilter[0]);
                        int VI_rawW = nb;
                        float *new_VI_raw = new float[savgolFilterW + VI_rawW - 1];
                        for (int m = 0; m < savgolFilterW + VI_rawW - 1; m++) {
                            if (m < (savgolFilterW - 1) / 2)
                                new_VI_raw[m] = VI_raw[0];
                            else if (m > ((savgolFilterW - 1) / 2 + VI_rawW - 1))
                                new_VI_raw[m] = VI_raw[VI_rawW - 1];
                            else
                                new_VI_raw[m] = VI_raw[m - (savgolFilterW - 1) / 2];
                        }
                        for (int m = 0; m < VI_rawW; m++) {
                            float temp = 0;
                            for (int n = 0; n < savgolFilterW; n++)
                                temp += savgolFilter[n] * new_VI_raw[m + n];
                            rst[m] = temp;
                        }
                        delete[] new_VI_raw;

                        // Calculate the weights for each point
                        float gdis = 0.0;
                        float *fl = new float[nb];
                        float maxdif = 0;
                        for (int k = 0; k < nb; k++) {
                            if ((VI_raw[k] - rst[k]) >= 0)
                                fl[k] = VI_raw[k] - rst[k];
                            else
                                fl[k] = rst[k] - VI_raw[k];

                            if (k == 0)
                                maxdif = fl[k];
                            else {
                                if (maxdif < fl[k])
                                    maxdif = fl[k];
                            }
                        }
                        for (int k = 0; k < nb; k++) {
                            if (VI_raw[k] >= rst[k]) {
                                fl[k] = 1.0;
                                gdis = gdis + fl[k] * (VI_raw[k] - rst[k]);
                            } else {
                                fl[k] = 1 - (rst[k] - VI_raw[k]) / maxdif;
                                gdis = gdis + fl[k] * (rst[k] - VI_raw[k]);
                            }
                        }

                        float *ra4 = new float[nb];
                        float *pre = new float[nb];
                        float ormax = gdis;
                        int loop_times = 1;
                        while (gdis <= ormax && loop_times < 15) {
                            loop_times += 1;
                            for (int k = 0; k < nb; k++) {
                                ra4[k] = (VI_raw[k] >= rst[k]) ? VI_raw[k] : rst[k];
                                pre[k] = rst[k];
                            }
                            // The Savitzky - Golay fitting
                            // savgolFilter = SAVGOL(4, 4, 0, 6); set the window width(4, 4) and degree(6) for repetition
                            double savgolFilter[] = {-0.00543880, 0.0435097, -0.152289, 0.304585, 0.619267, 0.304585, -0.152289, 0.0435097, -0.00543880};
                            int savgolFilterW = sizeof(savgolFilter) / sizeof(savgolFilter[0]);
                            int ra4W = nb;
                            float *new_ra4 = new float[savgolFilterW + ra4W - 1];
                            for (int m = 0; m < savgolFilterW + ra4W - 1; m++) {
                                if (m < (savgolFilterW - 1) / 2)
                                    new_ra4[m] = ra4[0];
                                else if (m > ((savgolFilterW - 1) / 2 + ra4W - 1))
                                    new_ra4[m] = ra4[ra4W - 1];
                                else
                                    new_ra4[m] = ra4[m - (savgolFilterW - 1) / 2];
                            }
                            for (int m = 0; m < ra4W; m++) {
                                float temp = 0;
                                for (int n = 0; n < savgolFilterW; n++)
                                    temp += savgolFilter[n] * new_ra4[m + n];
                                rst[m] = temp;
                            }
                            delete[] new_ra4;
                            ormax = gdis;
                            // Calculate the fitting - effect index
                            gdis = 0.0;
                            for (int k = 0; k < nb; k++) {
                                if ((VI_raw[k] - rst[k]) >= 0)
                                    gdis = gdis + fl[k] * (VI_raw[k] - rst[k]);
                                else
                                    gdis = gdis + fl[k] * (rst[k] - VI_raw[k]);
                            }
                        }
                        delete[] rst;
                        delete[] fl;
                        delete[] ra4;

                        for (int k = 0; k < nb; k++)
                            vector_out[i - win + (j - win) * (ns - 2 * win) + k * (ns - 2 * win) * (nl - 2 * win) + y * (ns - 2 * win) * (nl - 2 * win) * nb] = pre[k];
                        delete[] pre;
                    }
                } else {
                    for (int k = 0; k < nb; k++)
                        vector_out[i - win + (j - win) * (ns - 2 * win) + k * (ns - 2 * win) * (nl - 2 * win) + y * (ns - 2 * win) * (nl - 2 * win) * nb] = 0;
                }
                delete[] VI_raw;
                delete[] vector_QA;
                delete[] res_VI_raw;
            }
            delete[] similar_index;
            delete[] slope_intercept;
            delete[] corr_similar;
        }
    }
}
