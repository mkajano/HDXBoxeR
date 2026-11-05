# Returns standard deviation for uptake data for timecourses.

Calculates standard deviation for timecourse data.

## Usage

``` r
sd_timecourse(filepath)
```

## Arguments

- filepath:

  filepath to the All_results input file.

## Value

Data.frame with standard deviation.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
sd_timecourse(filepath=file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
#>    Protein.State Start End                Sequence Search.RT Charge
#> 1        Unbound    40  58     GPLGSKAVVPGPAEHPLQY      6.83      2
#> 2          bound    40  58     GPLGSKAVVPGPAEHPLQY      6.83      2
#> 3        Unbound    40  59    GPLGSKAVVPGPAEHPLQYN      6.44      2
#> 4          bound    40  59    GPLGSKAVVPGPAEHPLQYN      6.44      2
#> 5        Unbound    40  61  GPLGSKAVVPGPAEHPLQYNYT      7.10      2
#> 6        Unbound    40  62 GPLGSKAVVPGPAEHPLQYNYTF      8.54      2
#> 7        Unbound    47  63       VVPGPAEHPLQYNYTFW     10.86      2
#> 8          bound    47  63       VVPGPAEHPLQYNYTFW     10.86      2
#> 9        Unbound    62  77        FWYSRRTPGRPTSSQS      6.12      3
#> 10         bound    62  77        FWYSRRTPGRPTSSQS      6.12      3
#> 11       Unbound    63  77         WYSRRTPGRPTSSQS      5.33      3
#> 12         bound    63  77         WYSRRTPGRPTSSQS      5.33      3
#> 13       Unbound    64  77          YSRRTPGRPTSSQS      4.64      2
#> 14         bound    64  77          YSRRTPGRPTSSQS      4.64      2
#> 15       Unbound    78  88             YEQNIKQIGTF      7.61      2
#> 16         bound    78  88             YEQNIKQIGTF      7.61      2
#> 17       Unbound    82  88                 IKQIGTF      6.81      2
#> 18         bound    82  88                 IKQIGTF      6.81      2
#> 19       Unbound    89  97               ASVEQFWRF     10.72      2
#> 20         bound    89  97               ASVEQFWRF     10.72      2
#> 21       Unbound    95 107           WRFYSHMVRPGDL      8.44      3
#> 22         bound    95 107           WRFYSHMVRPGDL      8.44      3
#> 23       Unbound    95 113     WRFYSHMVRPGDLTGHSDF      8.44      3
#> 24       Unbound    95 115   WRFYSHMVRPGDLTGHSDFHL      8.65      3
#> 25         bound    95 115   WRFYSHMVRPGDLTGHSDFHL      8.65      3
#> 26       Unbound    98 107              YSHMVRPGDL      6.32      2
#> 27         bound    98 113        YSHMVRPGDLTGHSDF      6.72      2
#> 28       Unbound    98 115      YSHMVRPGDLTGHSDFHL      7.36      2
#> 29         bound    98 115      YSHMVRPGDLTGHSDFHL      7.36      2
#> 30       Unbound   114 124             HLFKEGIKPMW      8.07      2
#> 31         bound   114 124             HLFKEGIKPMW      8.07      2
#> 32       Unbound   114 135  HLFKEGIKPMWEDDANKNGGKW      7.63      3
#> 33         bound   114 135  HLFKEGIKPMWEDDANKNGGKW      7.63      3
#> 34       Unbound   116 124               FKEGIKPMW      7.98      2
#> 35         bound   116 124               FKEGIKPMW      7.98      2
#> 36       Unbound   116 135    FKEGIKPMWEDDANKNGGKW      7.38      3
#> 37         bound   116 135    FKEGIKPMWEDDANKNGGKW      7.38      3
#> 38       Unbound   136 147            IIRLRKGLASRC      5.79      3
#> 39         bound   136 147            IIRLRKGLASRC      5.79      3
#> 40       Unbound   136 149          IIRLRKGLASRCWE      6.83      3
#> 41         bound   136 149          IIRLRKGLASRCWE      6.83      3
#> 42       Unbound   140 147                RKGLASRC      4.88      2
#> 43         bound   140 147                RKGLASRC      4.88      2
#> 44       Unbound   154 161                AMLGEQFM     10.23      2
#> 45         bound   154 161                AMLGEQFM     10.23      2
#> 46       Unbound   169 177               AVVSVRFQE      6.71      2
#> 47         bound   169 177               AVVSVRFQE      6.71      2
#> 48       Unbound   176 183                QEDIISIW     11.81      2
#> 49         bound   176 183                QEDIISIW     11.81      2
#> 50       Unbound   176 193      QEDIISIWNKTASDQATT      9.39      2
#> 51         bound   176 193      QEDIISIWNKTASDQATT      9.39      2
#> 52       Unbound   178 193        DIISIWNKTASDQATT      9.28      2
#> 53         bound   178 193        DIISIWNKTASDQATT      9.28      2
#> 54       Unbound   179 193         IISIWNKTASDQATT      7.70      2
#> 55       Unbound   182 193            IWNKTASDQATT      5.51      2
#> 56         bound   182 193            IWNKTASDQATT      5.51      2
#> 57       Unbound   183 193             WNKTASDQATT      5.13      2
#> 58         bound   183 193             WNKTASDQATT      5.13      2
#> 59       Unbound   183 200      WNKTASDQATTARIRDTL      6.60      2
#> 60         bound   183 200      WNKTASDQATTARIRDTL      6.60      2
#> 61       Unbound   194 200                 ARIRDTL      5.75      2
#> 62         bound   194 200                 ARIRDTL      5.75      2
#> 63         bound   201 212            RRVLNLPPNTIM      7.56      2
#> 64       Unbound   201 213           RRVLNLPPNTIME      7.38      2
#> 65         bound   201 213           RRVLNLPPNTIME      7.38      2
#> 66       Unbound   213 220                EYKTHTDS      4.41      2
#> 67         bound   213 220                EYKTHTDS      4.41      2
#> 68       Unbound   213 225           EYKTHTDSIKMPG      5.65      2
#> 69         bound   213 225           EYKTHTDSIKMPG      5.65      2
#> 70       Unbound   214 220                 YKTHTDS      4.36      2
#> 71         bound   214 220                 YKTHTDS      4.36      2
#>    sd_t0.00s_X..Deut sd_t3.00s_X..Deut sd_t60.00s_X..Deut sd_t1800.00s_X..Deut
#> 1                 NA       0.138384729        0.070868423          0.041789153
#> 2                 NA       0.053799009        0.026950572          0.061760289
#> 3                 NA       0.154305973        0.072380937          0.062788003
#> 4                 NA       0.104423816        0.028884829          0.020550750
#> 5                 NA       0.159575478        0.083978172          0.054500765
#> 6                 NA       0.164658232        0.067099429          0.022869193
#> 7                 NA       0.108684559        0.036828431          0.023028967
#> 8                 NA       0.047648015        0.037819748          0.017578396
#> 9                 NA       0.130599898        0.054638204          0.056000000
#> 10                NA       0.081143905        0.046800997          0.043500958
#> 11                NA       0.220953238        0.087922314          0.049426713
#> 12                NA       0.013576941        0.245919363          0.036864617
#> 13                NA       0.056897569        0.051733290          0.098807894
#> 14                NA       0.037541089        0.032254199          0.056500737
#> 15                NA       0.022501852        0.016462078          0.030664855
#> 16                NA       0.047606022        0.013051181          0.023515952
#> 17                NA       0.071582121        0.039399662          0.021702534
#> 18                NA       0.059534304        0.039208843          0.034394767
#> 19                NA       0.041677332        0.013453624          0.014153916
#> 20                NA       0.061991935        0.025501634          0.036909800
#> 21                NA       0.053817593        0.016921387          0.069060360
#> 22                NA       0.084180362        0.028728615          0.028053520
#> 23                NA       0.117145209        0.123259888          0.011269428
#> 24                NA       0.117619443        0.087202829          0.016563011
#> 25                NA       0.078780285        0.045763887          0.070021425
#> 26                NA       0.002309401        0.020297783          0.026633312
#> 27                NA       0.031501323        0.019502137          0.035679126
#> 28                NA       0.100897638        0.056929781          0.072917762
#> 29                NA       0.085143017        0.071318534          0.028160256
#> 30                NA       0.057639685        0.019857828          0.028213472
#> 31                NA       0.067485801        0.038974351          0.023115651
#> 32                NA       0.033060551        0.189145447          0.009848858
#> 33                NA       0.314468334        0.158786020          0.108491167
#> 34                NA       0.029461840        0.009018500          0.028378396
#> 35                NA       0.026286879        0.017897858          0.051597804
#> 36                NA       0.082613558        0.091882171          0.036963946
#> 37                NA       0.115572488        0.088793769          0.064632293
#> 38                NA       0.073118625        0.062042996          0.064117080
#> 39                NA       0.114631293        0.050342163          0.066730802
#> 40                NA       0.024006943        0.055967252          0.064645185
#> 41                NA       0.131165290        0.013613719          0.038157568
#> 42                NA       0.105256829        0.035763109          0.043844422
#> 43                NA       0.090279197        0.072794231          0.014742230
#> 44                NA       0.018520259        0.003511885          0.024879711
#> 45                NA       0.102344191        0.035795717          0.011846237
#> 46                NA       0.030435725        0.020404248          0.003785939
#> 47                NA       0.019974984        0.021283797          0.025813433
#> 48                NA       0.017672955        0.028844410          0.024193663
#> 49                NA       0.006082763        0.007094599          0.071860977
#> 50                NA       0.074594459        0.053984566          0.117462051
#> 51                NA       0.062952363        0.066007575          0.153513300
#> 52                NA       0.013747727        0.032078030          0.078790862
#> 53                NA       0.075593210        0.061492547          0.034385074
#> 54                NA       0.074473709        0.076098620          0.043844422
#> 55                NA       0.000000000        0.000000000          0.000000000
#> 56                NA       0.000000000        0.000000000          0.000000000
#> 57                NA       0.074741778        0.073975221          0.062748174
#> 58                NA       0.060368314        0.057000000          0.001527525
#> 59                NA       0.078936683        0.272170902          0.055967252
#> 60                NA       0.262709599        0.056518434          0.302232251
#> 61                NA       0.032516662        0.048418316          0.030347982
#> 62                NA       0.054993939        0.007505553          0.046144700
#> 63                NA       0.075973680        0.047437678          0.031817186
#> 64                NA       0.074357246        0.046112182          0.046360903
#> 65                NA       0.076774562        0.027730849          0.022300972
#> 66                NA       0.046285347        0.042359572          0.070811958
#> 67                NA       0.124663547        0.102050641          0.023643181
#> 68                NA       0.099625967        0.028360771          0.042359572
#> 69                NA       0.019924859        0.027874720          0.020784610
#> 70                NA       0.026153394        0.024583192          0.014730920
#> 71                NA       0.024172988        0.014730920          0.050292478
#>    sd_t72000.00s_X..Deut sd_tFD_X..Deut
#> 1            0.033020196             NA
#> 2            0.060739883             NA
#> 3            0.056199051             NA
#> 4            0.069895160             NA
#> 5            0.086152965             NA
#> 6            0.030610456             NA
#> 7            0.022546249             NA
#> 8            0.017953644             NA
#> 9            0.040841156             NA
#> 10           0.065825527             NA
#> 11           0.037753587             NA
#> 12           0.074069787             NA
#> 13           0.079831072             NA
#> 14           0.001000000             NA
#> 15           0.033600595             NA
#> 16           0.038039453             NA
#> 17           0.035218366             NA
#> 18           0.027428695             NA
#> 19           0.014798649             NA
#> 20           0.030730007             NA
#> 21           0.040501029             NA
#> 22           0.048211340             NA
#> 23           0.078519637             NA
#> 24           0.027061658             NA
#> 25           0.045938365             NA
#> 26           0.030615900             NA
#> 27           0.046479386             NA
#> 28           0.037898989             NA
#> 29           0.046112182             NA
#> 30           0.047056703             NA
#> 31           0.027622455             NA
#> 32           0.126231269             NA
#> 33           0.077293812             NA
#> 34           0.037541089             NA
#> 35           0.063657940             NA
#> 36           0.082718398             NA
#> 37           0.005507571             NA
#> 38           0.027622455             NA
#> 39           0.045544850             NA
#> 40           0.025026652             NA
#> 41           0.099570745             NA
#> 42           0.078398554             NA
#> 43           0.059774577             NA
#> 44           0.011503623             NA
#> 45           0.053910420             NA
#> 46           0.037242449             NA
#> 47           0.025696952             NA
#> 48           0.024214321             NA
#> 49           0.031182260             NA
#> 50           0.050836339             NA
#> 51           0.078678671             NA
#> 52           0.055380502             NA
#> 53           0.072597061             NA
#> 54           0.054930259             NA
#> 55           0.000000000             NA
#> 56           0.000000000             NA
#> 57           0.066304852             NA
#> 58           0.073785726             NA
#> 59           0.104818891             NA
#> 60           0.030347982             NA
#> 61           0.032593455             NA
#> 62           0.019157244             NA
#> 63           0.020518285             NA
#> 64           0.045456939             NA
#> 65           0.056536124             NA
#> 66           0.056630381             NA
#> 67           0.093152205             NA
#> 68           0.030072135             NA
#> 69           0.029160476             NA
#> 70           0.023579652             NA
#> 71           0.072417769             NA
```
