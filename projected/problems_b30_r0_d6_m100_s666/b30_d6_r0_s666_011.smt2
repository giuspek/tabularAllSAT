(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
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
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A8))) (let ((.def_1 (and A2 .def_0))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A0))) (let ((.def_4 (not A5))) (let ((.def_5 (and .def_4 .def_3))) (let ((.def_6 (= .def_5 .def_2))) (let ((.def_7 (or A0 A13))) (let ((.def_8 (not .def_7))) (let ((.def_9 (not A12))) (let ((.def_10 (or A5 .def_9))) (let ((.def_11 (and .def_10 .def_8))) (let ((.def_12 (and .def_11 .def_6))) (let ((.def_13 (and A10 A14))) (let ((.def_14 (or .def_4 A0))) (let ((.def_15 (or .def_14 .def_13))) (let ((.def_16 (and .def_0 A3))) (let ((.def_17 (not .def_16))) (let ((.def_18 (and A16 A5))) (let ((.def_19 (and .def_18 .def_17))) (let ((.def_20 (or .def_19 .def_15))) (let ((.def_21 (or .def_20 .def_12))) (let ((.def_22 (not A3))) (let ((.def_23 (or .def_22 A7))) (let ((.def_24 (not .def_23))) (let ((.def_25 (not A4))) (let ((.def_26 (not A26))) (let ((.def_27 (or .def_26 .def_25))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and .def_28 .def_24))) (let ((.def_30 (not .def_29))) (let ((.def_31 (not A19))) (let ((.def_32 (and .def_0 .def_31))) (let ((.def_33 (not .def_32))) (let ((.def_34 (not A20))) (let ((.def_35 (not A25))) (let ((.def_36 (or .def_35 .def_34))) (let ((.def_37 (not .def_36))) (let ((.def_38 (and .def_37 .def_33))) (let ((.def_39 (and .def_38 .def_30))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A2))) (let ((.def_42 (and .def_35 .def_41))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A13))) (let ((.def_45 (and .def_44 A16))) (let ((.def_46 (or .def_45 .def_43))) (let ((.def_47 (not .def_46))) (let ((.def_48 (and A10 .def_22))) (let ((.def_49 (not .def_48))) (let ((.def_50 (not A23))) (let ((.def_51 (not A18))) (let ((.def_52 (and .def_51 .def_50))) (let ((.def_53 (and .def_52 .def_49))) (let ((.def_54 (not .def_53))) (let ((.def_55 (or .def_54 .def_47))) (let ((.def_56 (or .def_55 .def_40))) (let ((.def_57 (not .def_56))) (let ((.def_58 (= .def_57 .def_21))) (let ((.def_59 (not A1))) (let ((.def_60 (not A9))) (let ((.def_61 (or .def_60 .def_59))) (let ((.def_62 (not .def_61))) (let ((.def_63 (= .def_9 A4))) (let ((.def_64 (not .def_63))) (let ((.def_65 (or .def_64 .def_62))) (let ((.def_66 (and A28 A8))) (let ((.def_67 (not A6))) (let ((.def_68 (and A15 .def_67))) (let ((.def_69 (or .def_68 .def_66))) (let ((.def_70 (not .def_69))) (let ((.def_71 (or .def_70 .def_65))) (let ((.def_72 (not .def_71))) (let ((.def_73 (not A21))) (let ((.def_74 (or .def_60 .def_73))) (let ((.def_75 (not A29))) (let ((.def_76 (and .def_75 .def_59))) (let ((.def_77 (or .def_76 .def_74))) (let ((.def_78 (or .def_41 .def_44))) (let ((.def_79 (not .def_78))) (let ((.def_80 (or A22 A12))) (let ((.def_81 (not .def_80))) (let ((.def_82 (and .def_81 .def_79))) (let ((.def_83 (or .def_82 .def_77))) (let ((.def_84 (or .def_83 .def_72))) (let ((.def_85 (not A10))) (let ((.def_86 (or A3 .def_85))) (let ((.def_87 (not .def_86))) (let ((.def_88 (not A14))) (let ((.def_89 (or .def_35 .def_88))) (let ((.def_90 (or .def_89 .def_87))) (let ((.def_91 (not .def_90))) (let ((.def_92 (and .def_31 A23))) (let ((.def_93 (not .def_92))) (let ((.def_94 (and .def_51 A18))) (let ((.def_95 (or .def_94 .def_93))) (let ((.def_96 (not .def_95))) (let ((.def_97 (or .def_96 .def_91))) (let ((.def_98 (or A1 A19))) (let ((.def_99 (not .def_98))) (let ((.def_100 (not A27))) (let ((.def_101 (and A29 .def_100))) (let ((.def_102 (not .def_101))) (let ((.def_103 (and .def_102 .def_99))) (let ((.def_104 (= .def_50 A17))) (let ((.def_105 (not .def_104))) (let ((.def_106 (or A19 A6))) (let ((.def_107 (not .def_106))) (let ((.def_108 (= .def_107 .def_105))) (let ((.def_109 (not .def_108))) (let ((.def_110 (= .def_109 .def_103))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or .def_111 .def_97))) (let ((.def_113 (not .def_112))) (let ((.def_114 (and .def_113 .def_84))) (let ((.def_115 (not .def_114))) (let ((.def_116 (and .def_115 .def_58))) (let ((.def_117 (not .def_116))) .def_117)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
