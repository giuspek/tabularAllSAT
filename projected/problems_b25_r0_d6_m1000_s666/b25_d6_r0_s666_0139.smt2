(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
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
(assert (let ((.def_0 (not A0))) (let ((.def_1 (or A16 .def_0))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A22))) (let ((.def_4 (or .def_3 A11))) (let ((.def_5 (not .def_4))) (let ((.def_6 (and .def_5 .def_2))) (let ((.def_7 (not .def_6))) (let ((.def_8 (not A2))) (let ((.def_9 (not A9))) (let ((.def_10 (= .def_9 .def_8))) (let ((.def_11 (not A12))) (let ((.def_12 (not A14))) (let ((.def_13 (or .def_12 .def_11))) (let ((.def_14 (and .def_13 .def_10))) (let ((.def_15 (not .def_14))) (let ((.def_16 (and .def_15 .def_7))) (let ((.def_17 (not A6))) (let ((.def_18 (or A0 .def_17))) (let ((.def_19 (not A21))) (let ((.def_20 (not A4))) (let ((.def_21 (and .def_20 .def_19))) (let ((.def_22 (or .def_21 .def_18))) (let ((.def_23 (not .def_22))) (let ((.def_24 (not A18))) (let ((.def_25 (or A16 .def_24))) (let ((.def_26 (not .def_25))) (let ((.def_27 (not A5))) (let ((.def_28 (or .def_27 .def_0))) (let ((.def_29 (or .def_28 .def_26))) (let ((.def_30 (not .def_29))) (let ((.def_31 (and .def_30 .def_23))) (let ((.def_32 (or .def_31 .def_16))) (let ((.def_33 (and .def_24 A12))) (let ((.def_34 (not .def_33))) (let ((.def_35 (or A9 .def_20))) (let ((.def_36 (not .def_35))) (let ((.def_37 (= .def_36 .def_34))) (let ((.def_38 (not .def_37))) (let ((.def_39 (and A16 A5))) (let ((.def_40 (not .def_39))) (let ((.def_41 (and .def_17 A19))) (let ((.def_42 (or .def_41 .def_40))) (let ((.def_43 (not .def_42))) (let ((.def_44 (or .def_43 .def_38))) (let ((.def_45 (not .def_44))) (let ((.def_46 (and .def_11 A4))) (let ((.def_47 (and .def_19 .def_27))) (let ((.def_48 (not .def_47))) (let ((.def_49 (or .def_48 .def_46))) (let ((.def_50 (or .def_12 A4))) (let ((.def_51 (not A16))) (let ((.def_52 (and A18 .def_51))) (let ((.def_53 (not .def_52))) (let ((.def_54 (and .def_53 .def_50))) (let ((.def_55 (and .def_54 .def_49))) (let ((.def_56 (and .def_55 .def_45))) (let ((.def_57 (or .def_56 .def_32))) (let ((.def_58 (not A11))) (let ((.def_59 (or A4 .def_58))) (let ((.def_60 (not .def_59))) (let ((.def_61 (or .def_12 A10))) (let ((.def_62 (not .def_61))) (let ((.def_63 (and .def_62 .def_60))) (let ((.def_64 (not .def_63))) (let ((.def_65 (not A7))) (let ((.def_66 (not A19))) (let ((.def_67 (and .def_66 .def_65))) (let ((.def_68 (not .def_67))) (let ((.def_69 (and .def_58 A9))) (let ((.def_70 (= .def_69 .def_68))) (let ((.def_71 (or .def_70 .def_64))) (let ((.def_72 (not .def_71))) (let ((.def_73 (or A21 .def_12))) (let ((.def_74 (not .def_73))) (let ((.def_75 (and A3 A3))) (let ((.def_76 (not .def_75))) (let ((.def_77 (or .def_76 .def_74))) (let ((.def_78 (not .def_77))) (let ((.def_79 (and A17 A13))) (let ((.def_80 (not .def_79))) (let ((.def_81 (not A1))) (let ((.def_82 (and .def_81 A14))) (let ((.def_83 (not .def_82))) (let ((.def_84 (or .def_83 .def_80))) (let ((.def_85 (not .def_84))) (let ((.def_86 (= .def_85 .def_78))) (let ((.def_87 (not .def_86))) (let ((.def_88 (and .def_87 .def_72))) (let ((.def_89 (and .def_19 .def_0))) (let ((.def_90 (not .def_89))) (let ((.def_91 (not A20))) (let ((.def_92 (= .def_0 .def_91))) (let ((.def_93 (not .def_92))) (let ((.def_94 (or .def_93 .def_90))) (let ((.def_95 (not .def_94))) (let ((.def_96 (= .def_12 .def_0))) (let ((.def_97 (and A23 A11))) (let ((.def_98 (not .def_97))) (let ((.def_99 (and .def_98 .def_96))) (let ((.def_100 (or .def_99 .def_95))) (let ((.def_101 (and A17 A10))) (let ((.def_102 (not .def_101))) (let ((.def_103 (and A14 .def_24))) (let ((.def_104 (not .def_103))) (let ((.def_105 (or .def_104 .def_102))) (let ((.def_106 (or A21 A1))) (let ((.def_107 (not A15))) (let ((.def_108 (or .def_107 .def_27))) (let ((.def_109 (not .def_108))) (let ((.def_110 (and .def_109 .def_106))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or .def_111 .def_105))) (let ((.def_113 (not .def_112))) (let ((.def_114 (or .def_113 .def_100))) (let ((.def_115 (not .def_114))) (let ((.def_116 (or .def_115 .def_88))) (let ((.def_117 (and .def_116 .def_57))) (let ((.def_118 (not .def_117))) .def_118))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
