(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
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
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A5))) (let ((.def_1 (and A5 .def_0))) (let ((.def_2 (not .def_1))) (let ((.def_3 (not A20))) (let ((.def_4 (and .def_3 A20))) (let ((.def_5 (or .def_4 .def_2))) (let ((.def_6 (not .def_5))) (let ((.def_7 (or A1 A8))) (let ((.def_8 (not .def_7))) (let ((.def_9 (and A9 A10))) (let ((.def_10 (or .def_9 .def_8))) (let ((.def_11 (or .def_10 .def_6))) (let ((.def_12 (not A24))) (let ((.def_13 (or A22 .def_12))) (let ((.def_14 (= .def_3 A21))) (let ((.def_15 (not .def_14))) (let ((.def_16 (and .def_15 .def_13))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A23))) (let ((.def_19 (not A7))) (let ((.def_20 (= .def_19 .def_18))) (let ((.def_21 (and A8 A21))) (let ((.def_22 (not .def_21))) (let ((.def_23 (or .def_22 .def_20))) (let ((.def_24 (and .def_23 .def_17))) (let ((.def_25 (or .def_24 .def_11))) (let ((.def_26 (not .def_25))) (let ((.def_27 (and A22 A5))) (let ((.def_28 (not A16))) (let ((.def_29 (and .def_28 A11))) (let ((.def_30 (not .def_29))) (let ((.def_31 (or .def_30 .def_27))) (let ((.def_32 (not .def_31))) (let ((.def_33 (or A4 A23))) (let ((.def_34 (not .def_33))) (let ((.def_35 (not A18))) (let ((.def_36 (and A11 .def_35))) (let ((.def_37 (or .def_36 .def_34))) (let ((.def_38 (not .def_37))) (let ((.def_39 (or .def_38 .def_32))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A13))) (let ((.def_42 (and A18 .def_41))) (let ((.def_43 (or A18 .def_19))) (let ((.def_44 (and .def_43 .def_42))) (let ((.def_45 (not .def_44))) (let ((.def_46 (= .def_35 A7))) (let ((.def_47 (not .def_46))) (let ((.def_48 (not A6))) (let ((.def_49 (not A1))) (let ((.def_50 (or .def_49 .def_48))) (let ((.def_51 (= .def_50 .def_47))) (let ((.def_52 (and .def_51 .def_45))) (let ((.def_53 (not .def_52))) (let ((.def_54 (or .def_53 .def_40))) (let ((.def_55 (not .def_54))) (let ((.def_56 (and .def_55 .def_26))) (let ((.def_57 (or A4 A6))) (let ((.def_58 (not A11))) (let ((.def_59 (and .def_58 A16))) (let ((.def_60 (and .def_59 .def_57))) (let ((.def_61 (not A17))) (let ((.def_62 (not A0))) (let ((.def_63 (or .def_62 .def_61))) (let ((.def_64 (not .def_63))) (let ((.def_65 (and A16 A21))) (let ((.def_66 (and .def_65 .def_64))) (let ((.def_67 (not .def_66))) (let ((.def_68 (and .def_67 .def_60))) (let ((.def_69 (not .def_68))) (let ((.def_70 (not A9))) (let ((.def_71 (or A4 .def_70))) (let ((.def_72 (not .def_71))) (let ((.def_73 (not A21))) (let ((.def_74 (not A12))) (let ((.def_75 (= .def_74 .def_73))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and .def_76 .def_72))) (let ((.def_78 (not .def_77))) (let ((.def_79 (or A21 A11))) (let ((.def_80 (or .def_58 A16))) (let ((.def_81 (and .def_80 .def_79))) (let ((.def_82 (or .def_81 .def_78))) (let ((.def_83 (not .def_82))) (let ((.def_84 (and .def_83 .def_69))) (let ((.def_85 (and A17 A23))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and A1 .def_62))) (let ((.def_88 (not .def_87))) (let ((.def_89 (or .def_88 .def_86))) (let ((.def_90 (not .def_89))) (let ((.def_91 (or .def_19 A10))) (let ((.def_92 (not A10))) (let ((.def_93 (and .def_92 A7))) (let ((.def_94 (not .def_93))) (let ((.def_95 (and .def_94 .def_91))) (let ((.def_96 (or .def_95 .def_90))) (let ((.def_97 (or A13 .def_41))) (let ((.def_98 (not .def_97))) (let ((.def_99 (and A21 .def_62))) (let ((.def_100 (or .def_99 .def_98))) (let ((.def_101 (or .def_41 .def_0))) (let ((.def_102 (and .def_61 A8))) (let ((.def_103 (and .def_102 .def_101))) (let ((.def_104 (or .def_103 .def_100))) (let ((.def_105 (not .def_104))) (let ((.def_106 (or .def_105 .def_96))) (let ((.def_107 (not .def_106))) (let ((.def_108 (and .def_107 .def_84))) (let ((.def_109 (and .def_108 .def_56))) (let ((.def_110 (not .def_109))) .def_110))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)