(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (or A26 A8))) (let ((.def_1 (not A15))) (let ((.def_2 (= .def_1 A13))) (let ((.def_3 (= .def_2 .def_0))) (let ((.def_4 (or A21 A8))) (let ((.def_5 (not A9))) (let ((.def_6 (or A1 .def_5))) (let ((.def_7 (and .def_6 .def_4))) (let ((.def_8 (not .def_7))) (let ((.def_9 (and .def_8 .def_3))) (let ((.def_10 (not A21))) (let ((.def_11 (not A5))) (let ((.def_12 (or .def_11 .def_10))) (let ((.def_13 (not .def_12))) (let ((.def_14 (or .def_10 A21))) (let ((.def_15 (= .def_14 .def_13))) (let ((.def_16 (not A7))) (let ((.def_17 (not A13))) (let ((.def_18 (and .def_17 .def_16))) (let ((.def_19 (not A25))) (let ((.def_20 (or A19 .def_19))) (let ((.def_21 (and .def_20 .def_18))) (let ((.def_22 (and .def_21 .def_15))) (let ((.def_23 (and .def_22 .def_9))) (let ((.def_24 (not .def_23))) (let ((.def_25 (not A24))) (let ((.def_26 (or .def_11 .def_25))) (let ((.def_27 (not A29))) (let ((.def_28 (or A7 .def_27))) (let ((.def_29 (not .def_28))) (let ((.def_30 (and .def_29 .def_26))) (let ((.def_31 (not .def_30))) (let ((.def_32 (not A23))) (let ((.def_33 (or .def_32 A22))) (let ((.def_34 (not .def_33))) (let ((.def_35 (not A19))) (let ((.def_36 (= .def_35 A5))) (let ((.def_37 (and .def_36 .def_34))) (let ((.def_38 (not .def_37))) (let ((.def_39 (and .def_38 .def_31))) (let ((.def_40 (not .def_39))) (let ((.def_41 (not A28))) (let ((.def_42 (or .def_11 .def_41))) (let ((.def_43 (and A8 A14))) (let ((.def_44 (not .def_43))) (let ((.def_45 (and .def_44 .def_42))) (let ((.def_46 (not A6))) (let ((.def_47 (and .def_25 .def_46))) (let ((.def_48 (not A22))) (let ((.def_49 (and .def_48 A23))) (let ((.def_50 (and .def_49 .def_47))) (let ((.def_51 (not .def_50))) (let ((.def_52 (or .def_51 .def_45))) (let ((.def_53 (not .def_52))) (let ((.def_54 (or .def_53 .def_40))) (let ((.def_55 (and .def_54 .def_24))) (let ((.def_56 (not .def_55))) (let ((.def_57 (not A0))) (let ((.def_58 (or .def_41 .def_57))) (let ((.def_59 (not A14))) (let ((.def_60 (not A4))) (let ((.def_61 (and .def_60 .def_59))) (let ((.def_62 (and .def_61 .def_58))) (let ((.def_63 (and .def_25 .def_19))) (let ((.def_64 (not .def_63))) (let ((.def_65 (and A15 .def_19))) (let ((.def_66 (or .def_65 .def_64))) (let ((.def_67 (and .def_66 .def_62))) (let ((.def_68 (and A9 .def_59))) (let ((.def_69 (not .def_68))) (let ((.def_70 (and .def_17 A16))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_71 .def_69))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and A28 .def_11))) (let ((.def_75 (and .def_5 .def_48))) (let ((.def_76 (not .def_75))) (let ((.def_77 (or .def_76 .def_74))) (let ((.def_78 (not .def_77))) (let ((.def_79 (and .def_78 .def_73))) (let ((.def_80 (not .def_79))) (let ((.def_81 (and .def_80 .def_67))) (let ((.def_82 (not .def_81))) (let ((.def_83 (not A18))) (let ((.def_84 (or A27 .def_83))) (let ((.def_85 (not .def_84))) (let ((.def_86 (or A28 A20))) (let ((.def_87 (and .def_86 .def_85))) (let ((.def_88 (not .def_87))) (let ((.def_89 (not A12))) (let ((.def_90 (= .def_89 .def_83))) (let ((.def_91 (not A2))) (let ((.def_92 (or A25 .def_91))) (let ((.def_93 (or .def_92 .def_90))) (let ((.def_94 (and .def_93 .def_88))) (let ((.def_95 (and A23 .def_41))) (let ((.def_96 (not .def_95))) (let ((.def_97 (= A16 A23))) (let ((.def_98 (or .def_97 .def_96))) (let ((.def_99 (not .def_98))) (let ((.def_100 (and A23 .def_11))) (let ((.def_101 (not .def_100))) (let ((.def_102 (or .def_91 .def_46))) (let ((.def_103 (and .def_102 .def_101))) (let ((.def_104 (or .def_103 .def_99))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_94))) (let ((.def_107 (= .def_106 .def_82))) (let ((.def_108 (and .def_107 .def_56))) (let ((.def_109 (not .def_108))) .def_109)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)