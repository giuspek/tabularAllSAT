(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(assert (let ((.def_0 (not A15))) (let ((.def_1 (not A7))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not A16))) (let ((.def_4 (not A10))) (let ((.def_5 (or .def_4 .def_3))) (let ((.def_6 (not .def_5))) (let ((.def_7 (and .def_6 .def_2))) (let ((.def_8 (not .def_7))) (let ((.def_9 (or A8 A7))) (let ((.def_10 (not A4))) (let ((.def_11 (not A3))) (let ((.def_12 (and .def_11 .def_10))) (let ((.def_13 (not .def_12))) (let ((.def_14 (or .def_13 .def_9))) (let ((.def_15 (not .def_14))) (let ((.def_16 (and .def_15 .def_8))) (let ((.def_17 (not .def_16))) (let ((.def_18 (and A1 A20))) (let ((.def_19 (not .def_18))) (let ((.def_20 (not A8))) (let ((.def_21 (not A9))) (let ((.def_22 (= .def_21 .def_20))) (let ((.def_23 (not .def_22))) (let ((.def_24 (= .def_23 .def_19))) (let ((.def_25 (not .def_24))) (let ((.def_26 (or A17 A14))) (let ((.def_27 (not .def_26))) (let ((.def_28 (not A21))) (let ((.def_29 (or .def_28 A1))) (let ((.def_30 (not .def_29))) (let ((.def_31 (and .def_30 .def_27))) (let ((.def_32 (and .def_31 .def_25))) (let ((.def_33 (not .def_32))) (let ((.def_34 (and .def_33 .def_17))) (let ((.def_35 (not .def_34))) (let ((.def_36 (not A0))) (let ((.def_37 (= A12 .def_36))) (let ((.def_38 (or A12 A17))) (let ((.def_39 (not .def_38))) (let ((.def_40 (and .def_39 .def_37))) (let ((.def_41 (not .def_40))) (let ((.def_42 (and A8 A15))) (let ((.def_43 (= A23 A5))) (let ((.def_44 (or .def_43 .def_42))) (let ((.def_45 (not .def_44))) (let ((.def_46 (or .def_45 .def_41))) (let ((.def_47 (not .def_46))) (let ((.def_48 (or .def_3 A21))) (let ((.def_49 (or .def_20 A13))) (let ((.def_50 (not .def_49))) (let ((.def_51 (or .def_50 .def_48))) (let ((.def_52 (not .def_51))) (let ((.def_53 (and A19 A22))) (let ((.def_54 (and .def_28 A20))) (let ((.def_55 (not .def_54))) (let ((.def_56 (and .def_55 .def_53))) (let ((.def_57 (or .def_56 .def_52))) (let ((.def_58 (not .def_57))) (let ((.def_59 (or .def_58 .def_47))) (let ((.def_60 (not .def_59))) (let ((.def_61 (or .def_60 .def_35))) (let ((.def_62 (not .def_61))) (let ((.def_63 (not A22))) (let ((.def_64 (not A14))) (let ((.def_65 (and .def_64 .def_63))) (let ((.def_66 (not .def_65))) (let ((.def_67 (or A14 A14))) (let ((.def_68 (not .def_67))) (let ((.def_69 (and .def_68 .def_66))) (let ((.def_70 (not A17))) (let ((.def_71 (not A12))) (let ((.def_72 (= .def_71 .def_70))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_3 .def_28))) (let ((.def_75 (not .def_74))) (let ((.def_76 (and .def_75 .def_73))) (let ((.def_77 (and .def_76 .def_69))) (let ((.def_78 (not .def_77))) (let ((.def_79 (not A18))) (let ((.def_80 (= .def_79 .def_4))) (let ((.def_81 (not .def_80))) (let ((.def_82 (and .def_28 .def_71))) (let ((.def_83 (not .def_82))) (let ((.def_84 (or .def_83 .def_81))) (let ((.def_85 (not .def_84))) (let ((.def_86 (and .def_21 A4))) (let ((.def_87 (and A14 A8))) (let ((.def_88 (and .def_87 .def_86))) (let ((.def_89 (not .def_88))) (let ((.def_90 (and .def_89 .def_85))) (let ((.def_91 (not .def_90))) (let ((.def_92 (or .def_91 .def_78))) (let ((.def_93 (or A3 .def_10))) (let ((.def_94 (not .def_93))) (let ((.def_95 (or .def_36 A20))) (let ((.def_96 (and .def_95 .def_94))) (let ((.def_97 (not .def_96))) (let ((.def_98 (and A12 .def_1))) (let ((.def_99 (not .def_98))) (let ((.def_100 (not A6))) (let ((.def_101 (and .def_79 .def_100))) (let ((.def_102 (and .def_101 .def_99))) (let ((.def_103 (and .def_102 .def_97))) (let ((.def_104 (not .def_103))) (let ((.def_105 (and .def_79 A17))) (let ((.def_106 (not .def_105))) (let ((.def_107 (= A0 .def_21))) (let ((.def_108 (or .def_107 .def_106))) (let ((.def_109 (not .def_108))) (let ((.def_110 (or .def_4 A18))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or A21 A21))) (let ((.def_113 (not .def_112))) (let ((.def_114 (or .def_113 .def_111))) (let ((.def_115 (not .def_114))) (let ((.def_116 (or .def_115 .def_109))) (let ((.def_117 (or .def_116 .def_104))) (let ((.def_118 (not .def_117))) (let ((.def_119 (and .def_118 .def_92))) (let ((.def_120 (not .def_119))) (let ((.def_121 (or .def_120 .def_62))) (let ((.def_122 (not .def_121))) .def_122))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
