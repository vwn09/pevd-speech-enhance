# Introduction
This is the code repository for the PEVD-based speech enhancement algorithm for [1] which uses the SMD [2]. The algorithm is multi-channel, blind and unsupervised, and does not require any prior information such as array geometry, source direction etc. It does not rely on any noise and relative transfer function estimates.

Demo page: https://vwn09.github.io/pevd-enhance/

Contact Vincent Neo if you have any questions: <vwn09@ic.ac.uk>

## Instructions
1. Clone the repository `git clone https://github.com/vwn09/pevd-speech-enhance.git`
2. Run `example`

## Code
`example` demonstrates on how `pevd_enhance` can be used. `pevd_enhance` is designed to be a self-contained 'lightweight' function that executes [1]. The minimum input for the function is a matrix of multi-channel audio (N samples, M channels) and the algorithm will generate decorrelated outputs (N' samples, M channels). The enhanced signal is in the first channel.

If you use this code in your research, please cite our work [1] and acknowledge [2].

## References
[1] V. W. Neo, C. Evers, and P. A. Naylor, "Enhancement of noisy reverberant speech using polynomial matrix eigenvalue decomposition," *IEEE/ACM Trans. Audio, Speech and Lang. Process.*, vol. 28, 2021, doi:10.1109/TASLP.2021.3120630

[2] S. Redif, S. Weiss, and J. G. McWhirter, "Sequential matrix diagonalisation algorithms for polynomial EVD of parahermitian matrices," *IEEE Trans. Signal Process.*, vol. 63, no. 1, pp. 81-89, Jan. 2015.
