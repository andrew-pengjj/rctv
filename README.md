# The code of RCTV

## RCTV model is:
$\min_{\mathcal{U},\mathbf{V},\mathcal{S},\mathcal{E}} \quad \Vert \mathcal{U} \Vert_*+\lambda \Vert \mathcal{S} \Vert_1+\beta \Vert \mathcal{E} \Vert_F^2$ 

 s.t. $\mathcal{Y} = \mathcal{U}\times_3 \mathbf{V} + \mathcal{S} + \mathcal{E}, \mathbf{V}^T\mathbf{V}=\mathbf{I}$.
### Run ``` demo.m ``` to show the denoised results on simu_indian with size of $145\times 145\times 220$.
| Metric      | Value |
|---------| --------- |
| MPSNR   | 46.81     |
| MSSIM   | 0.9938    |
| ERGAS   | 12.23     |
| Times   | 5.34(s)   |

If you use this code, please cite this [Paper](https://ieeexplore.ieee.org/abstract/document/9989343 "悬停显示")
