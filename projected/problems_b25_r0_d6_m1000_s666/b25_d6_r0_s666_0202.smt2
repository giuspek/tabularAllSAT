(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
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
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A14))) (let ((.def_1 (or .def_0 A0))) (let ((.def_2 (or A6 A23))) (let ((.def_3 (and .def_2 .def_1))) (let ((.def_4 (not A17))) (let ((.def_5 (or .def_4 A23))) (let ((.def_6 (not .def_5))) (let ((.def_7 (not A4))) (let ((.def_8 (= A16 .def_7))) (let ((.def_9 (= .def_8 .def_6))) (let ((.def_10 (not .def_9))) (let ((.def_11 (or .def_10 .def_3))) (let ((.def_12 (not A5))) (let ((.def_13 (or A11 .def_12))) (let ((.def_14 (not .def_13))) (let ((.def_15 (not A24))) (let ((.def_16 (and .def_7 .def_15))) (let ((.def_17 (not .def_16))) (let ((.def_18 (or .def_17 .def_14))) (let ((.def_19 (not .def_18))) (let ((.def_20 (or A18 A11))) (let ((.def_21 (not .def_20))) (let ((.def_22 (not A15))) (let ((.def_23 (not A16))) (let ((.def_24 (or .def_23 .def_22))) (let ((.def_25 (or .def_24 .def_21))) (let ((.def_26 (or .def_25 .def_19))) (let ((.def_27 (not .def_26))) (let ((.def_28 (or .def_27 .def_11))) (let ((.def_29 (not .def_28))) (let ((.def_30 (and A18 A15))) (let ((.def_31 (and A14 A24))) (let ((.def_32 (not .def_31))) (let ((.def_33 (and .def_32 .def_30))) (let ((.def_34 (not .def_33))) (let ((.def_35 (not A13))) (let ((.def_36 (and A22 .def_35))) (let ((.def_37 (not .def_36))) (let ((.def_38 (or A24 .def_15))) (let ((.def_39 (not .def_38))) (let ((.def_40 (or .def_39 .def_37))) (let ((.def_41 (not .def_40))) (let ((.def_42 (and .def_41 .def_34))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A0))) (let ((.def_45 (or .def_44 A0))) (let ((.def_46 (not .def_45))) (let ((.def_47 (not A11))) (let ((.def_48 (and A14 .def_47))) (let ((.def_49 (or .def_48 .def_46))) (let ((.def_50 (not .def_49))) (let ((.def_51 (or A4 .def_22))) (let ((.def_52 (not .def_51))) (let ((.def_53 (not A19))) (let ((.def_54 (or .def_53 A21))) (let ((.def_55 (not .def_54))) (let ((.def_56 (= .def_55 .def_52))) (let ((.def_57 (not .def_56))) (let ((.def_58 (or .def_57 .def_50))) (let ((.def_59 (and .def_58 .def_43))) (let ((.def_60 (not .def_59))) (let ((.def_61 (and .def_60 .def_29))) (let ((.def_62 (and A0 A2))) (let ((.def_63 (not A22))) (let ((.def_64 (and .def_63 .def_44))) (let ((.def_65 (or .def_64 .def_62))) (let ((.def_66 (not .def_65))) (let ((.def_67 (not A20))) (let ((.def_68 (and .def_67 A23))) (let ((.def_69 (not .def_68))) (let ((.def_70 (= A4 .def_4))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_71 .def_69))) (let ((.def_73 (and .def_72 .def_66))) (let ((.def_74 (not A7))) (let ((.def_75 (or A21 .def_74))) (let ((.def_76 (and A11 .def_67))) (let ((.def_77 (not .def_76))) (let ((.def_78 (or .def_77 .def_75))) (let ((.def_79 (and A17 A10))) (let ((.def_80 (not .def_79))) (let ((.def_81 (not A2))) (let ((.def_82 (and A14 .def_81))) (let ((.def_83 (or .def_82 .def_80))) (let ((.def_84 (not .def_83))) (let ((.def_85 (or .def_84 .def_78))) (let ((.def_86 (and .def_85 .def_73))) (let ((.def_87 (and A23 A17))) (let ((.def_88 (= .def_7 .def_35))) (let ((.def_89 (or .def_88 .def_87))) (let ((.def_90 (not A1))) (let ((.def_91 (not A6))) (let ((.def_92 (or .def_91 .def_90))) (let ((.def_93 (= .def_67 A4))) (let ((.def_94 (and .def_93 .def_92))) (let ((.def_95 (not .def_94))) (let ((.def_96 (and .def_95 .def_89))) (let ((.def_97 (not .def_96))) (let ((.def_98 (not A21))) (let ((.def_99 (and A18 .def_98))) (let ((.def_100 (not .def_99))) (let ((.def_101 (and .def_44 .def_63))) (let ((.def_102 (not .def_101))) (let ((.def_103 (and .def_102 .def_100))) (let ((.def_104 (not .def_103))) (let ((.def_105 (= A11 .def_12))) (let ((.def_106 (= .def_91 .def_98))) (let ((.def_107 (not .def_106))) (let ((.def_108 (= .def_107 .def_105))) (let ((.def_109 (or .def_108 .def_104))) (let ((.def_110 (= .def_109 .def_97))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or .def_111 .def_86))) (let ((.def_113 (not .def_112))) (let ((.def_114 (or .def_113 .def_61))) .def_114))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
