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
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(assert (let ((.def_0 (not A3))) (let ((.def_1 (not A21))) (let ((.def_2 (and .def_1 .def_0))) (let ((.def_3 (not .def_2))) (let ((.def_4 (not A14))) (let ((.def_5 (or .def_4 A14))) (let ((.def_6 (not .def_5))) (let ((.def_7 (or .def_6 .def_3))) (let ((.def_8 (not A5))) (let ((.def_9 (or A10 .def_8))) (let ((.def_10 (not A13))) (let ((.def_11 (and .def_10 A6))) (let ((.def_12 (and .def_11 .def_9))) (let ((.def_13 (or .def_12 .def_7))) (let ((.def_14 (not .def_13))) (let ((.def_15 (not A18))) (let ((.def_16 (and A16 .def_15))) (let ((.def_17 (and .def_15 A8))) (let ((.def_18 (not .def_17))) (let ((.def_19 (or .def_18 .def_16))) (let ((.def_20 (and A14 A0))) (let ((.def_21 (not A11))) (let ((.def_22 (and A9 .def_21))) (let ((.def_23 (= .def_22 .def_20))) (let ((.def_24 (not .def_23))) (let ((.def_25 (or .def_24 .def_19))) (let ((.def_26 (or .def_25 .def_14))) (let ((.def_27 (not .def_26))) (let ((.def_28 (and A21 A22))) (let ((.def_29 (not .def_28))) (let ((.def_30 (not A4))) (let ((.def_31 (and .def_30 A21))) (let ((.def_32 (and .def_31 .def_29))) (let ((.def_33 (not A10))) (let ((.def_34 (and .def_33 .def_4))) (let ((.def_35 (not .def_34))) (let ((.def_36 (and .def_4 .def_8))) (let ((.def_37 (and .def_36 .def_35))) (let ((.def_38 (not .def_37))) (let ((.def_39 (and .def_38 .def_32))) (let ((.def_40 (not A7))) (let ((.def_41 (not A0))) (let ((.def_42 (or .def_41 .def_40))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A12))) (let ((.def_45 (not A23))) (let ((.def_46 (or .def_45 .def_44))) (let ((.def_47 (and .def_46 .def_43))) (let ((.def_48 (not .def_47))) (let ((.def_49 (= .def_21 A22))) (let ((.def_50 (not .def_49))) (let ((.def_51 (not A9))) (let ((.def_52 (not A1))) (let ((.def_53 (or .def_52 .def_51))) (let ((.def_54 (not .def_53))) (let ((.def_55 (or .def_54 .def_50))) (let ((.def_56 (not .def_55))) (let ((.def_57 (= .def_56 .def_48))) (let ((.def_58 (not .def_57))) (let ((.def_59 (and .def_58 .def_39))) (let ((.def_60 (or .def_59 .def_27))) (let ((.def_61 (not .def_60))) (let ((.def_62 (and .def_52 .def_40))) (let ((.def_63 (not .def_62))) (let ((.def_64 (not A2))) (let ((.def_65 (or .def_64 .def_64))) (let ((.def_66 (not .def_65))) (let ((.def_67 (and .def_66 .def_63))) (let ((.def_68 (not .def_67))) (let ((.def_69 (and .def_45 .def_52))) (let ((.def_70 (or A16 A2))) (let ((.def_71 (not .def_70))) (let ((.def_72 (= .def_71 .def_69))) (let ((.def_73 (not .def_72))) (let ((.def_74 (or .def_73 .def_68))) (let ((.def_75 (or .def_30 A18))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and A17 A14))) (let ((.def_78 (and .def_77 .def_76))) (let ((.def_79 (not .def_78))) (let ((.def_80 (or .def_21 A23))) (let ((.def_81 (or A17 A1))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and .def_82 .def_80))) (let ((.def_84 (or .def_83 .def_79))) (let ((.def_85 (or .def_84 .def_74))) (let ((.def_86 (or A21 .def_44))) (let ((.def_87 (or .def_41 A11))) (let ((.def_88 (and .def_87 .def_86))) (let ((.def_89 (or .def_4 A0))) (let ((.def_90 (not .def_89))) (let ((.def_91 (and .def_8 .def_64))) (let ((.def_92 (not .def_91))) (let ((.def_93 (and .def_92 .def_90))) (let ((.def_94 (or .def_93 .def_88))) (let ((.def_95 (not .def_94))) (let ((.def_96 (and .def_0 A22))) (let ((.def_97 (or A11 A5))) (let ((.def_98 (not .def_97))) (let ((.def_99 (or .def_98 .def_96))) (let ((.def_100 (and A23 .def_51))) (let ((.def_101 (or A3 A11))) (let ((.def_102 (not .def_101))) (let ((.def_103 (= .def_102 .def_100))) (let ((.def_104 (not .def_103))) (let ((.def_105 (or .def_104 .def_99))) (let ((.def_106 (or .def_105 .def_95))) (let ((.def_107 (= .def_106 .def_85))) (let ((.def_108 (and .def_107 .def_61))) .def_108))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
