(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A19 () Bool)
(declare-fun A23 () Bool)
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A7))) (let ((.def_1 (or A3 .def_0))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A16))) (let ((.def_4 (or .def_3 A15))) (let ((.def_5 (or .def_4 .def_2))) (let ((.def_6 (not .def_5))) (let ((.def_7 (or A10 A6))) (let ((.def_8 (and A4 A15))) (let ((.def_9 (or .def_8 .def_7))) (let ((.def_10 (not .def_9))) (let ((.def_11 (and .def_10 .def_6))) (let ((.def_12 (not .def_11))) (let ((.def_13 (not A8))) (let ((.def_14 (= .def_13 A9))) (let ((.def_15 (not .def_14))) (let ((.def_16 (not A27))) (let ((.def_17 (and .def_16 A13))) (let ((.def_18 (and .def_17 .def_15))) (let ((.def_19 (not .def_18))) (let ((.def_20 (not A19))) (let ((.def_21 (not A3))) (let ((.def_22 (and .def_21 .def_20))) (let ((.def_23 (not A12))) (let ((.def_24 (not A14))) (let ((.def_25 (and .def_24 .def_23))) (let ((.def_26 (not .def_25))) (let ((.def_27 (and .def_26 .def_22))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and .def_28 .def_19))) (let ((.def_30 (not .def_29))) (let ((.def_31 (or .def_30 .def_12))) (let ((.def_32 (not A13))) (let ((.def_33 (= A0 .def_32))) (let ((.def_34 (not .def_33))) (let ((.def_35 (not A29))) (let ((.def_36 (and A7 .def_35))) (let ((.def_37 (not .def_36))) (let ((.def_38 (and .def_37 .def_34))) (let ((.def_39 (not A11))) (let ((.def_40 (and A19 .def_39))) (let ((.def_41 (not .def_40))) (let ((.def_42 (or A29 .def_21))) (let ((.def_43 (not .def_42))) (let ((.def_44 (= .def_43 .def_41))) (let ((.def_45 (not .def_44))) (let ((.def_46 (and .def_45 .def_38))) (let ((.def_47 (= A14 A10))) (let ((.def_48 (not .def_47))) (let ((.def_49 (or A25 A23))) (let ((.def_50 (not .def_49))) (let ((.def_51 (or .def_50 .def_48))) (let ((.def_52 (not A26))) (let ((.def_53 (not A17))) (let ((.def_54 (and .def_53 .def_52))) (let ((.def_55 (and .def_35 A19))) (let ((.def_56 (and .def_55 .def_54))) (let ((.def_57 (not .def_56))) (let ((.def_58 (or .def_57 .def_51))) (let ((.def_59 (and .def_58 .def_46))) (let ((.def_60 (or .def_59 .def_31))) (let ((.def_61 (not A15))) (let ((.def_62 (and .def_20 .def_61))) (let ((.def_63 (not A6))) (let ((.def_64 (or .def_16 .def_63))) (let ((.def_65 (and .def_64 .def_62))) (let ((.def_66 (not .def_65))) (let ((.def_67 (and .def_23 .def_24))) (let ((.def_68 (and A6 .def_13))) (let ((.def_69 (or .def_68 .def_67))) (let ((.def_70 (not .def_69))) (let ((.def_71 (and .def_70 .def_66))) (let ((.def_72 (and A6 A17))) (let ((.def_73 (not .def_72))) (let ((.def_74 (not A4))) (let ((.def_75 (or A8 .def_74))) (let ((.def_76 (not .def_75))) (let ((.def_77 (or .def_76 .def_73))) (let ((.def_78 (not A2))) (let ((.def_79 (and .def_78 A9))) (let ((.def_80 (or .def_3 A29))) (let ((.def_81 (not .def_80))) (let ((.def_82 (and .def_81 .def_79))) (let ((.def_83 (or .def_82 .def_77))) (let ((.def_84 (= .def_83 .def_71))) (let ((.def_85 (or .def_24 A29))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and .def_24 .def_78))) (let ((.def_88 (not .def_87))) (let ((.def_89 (= .def_88 .def_86))) (let ((.def_90 (not A25))) (let ((.def_91 (or .def_78 .def_90))) (let ((.def_92 (not A5))) (let ((.def_93 (and .def_92 A13))) (let ((.def_94 (or .def_93 .def_91))) (let ((.def_95 (not .def_94))) (let ((.def_96 (or .def_95 .def_89))) (let ((.def_97 (not .def_96))) (let ((.def_98 (and A12 .def_0))) (let ((.def_99 (or .def_23 A29))) (let ((.def_100 (and .def_99 .def_98))) (let ((.def_101 (not .def_100))) (let ((.def_102 (or .def_53 .def_16))) (let ((.def_103 (not .def_102))) (let ((.def_104 (= A9 A9))) (let ((.def_105 (and .def_104 .def_103))) (let ((.def_106 (not .def_105))) (let ((.def_107 (and .def_106 .def_101))) (let ((.def_108 (or .def_107 .def_97))) (let ((.def_109 (or .def_108 .def_84))) (let ((.def_110 (not .def_109))) (let ((.def_111 (or .def_110 .def_60))) .def_111)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
