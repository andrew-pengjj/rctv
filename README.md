# The code of RCTV

>> ## RCTV model is:
>>>> $\min_{\mathcal{U},\mathbf{V},\mathcal{S},\mathcal{E}} \quad \Vert \mathcal{U} \Vert_*+\lambda \Vert \mathcal{S} \Vert_1+\beta \Vert \mathcal{E} \Vert_F^2$ 

>>>>>> s.t. $\mathcal{Y} = \mathcal{U}\times_3 \mathbf{V} + \mathcal{S} + \mathcal{E}, \mathbf{V}^T\mathbf{V}=\mathbf{I}$.
>> ### Run ``` demo.m ``` to show the denoised results on simu_indian with 
>> ### The denoised resultsis:
>>    MPSNR: 46.81, MSSIM:0.9938, ERGAS: 12.23, Times: 5.34s 

| Metric      | Value |
|---------| --------- |
| MPSNR   | 46.81     |
| MSSIM   | 0.9938    |
| ERGAS   | 12.23     |
| Times   | 5.34(s)   |

>>> For more information, please see the [Paper](https://ieeexplore.ieee.org/abstract/document/9989343 "悬停显示")
