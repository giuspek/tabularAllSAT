(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A23))) (let ((.def_1 (not A16))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not .def_2))) (let ((.def_4 (not A4))) (let ((.def_5 (or A23 .def_4))) (let ((.def_6 (or .def_5 .def_3))) (let ((.def_7 (or A3 A10))) (let ((.def_8 (or A1 A16))) (let ((.def_9 (and .def_8 .def_7))) (let ((.def_10 (not .def_9))) (let ((.def_11 (and .def_10 .def_6))) (let ((.def_12 (not .def_11))) (let ((.def_13 (not A5))) (let ((.def_14 (or A22 .def_13))) (let ((.def_15 (not .def_14))) (let ((.def_16 (not A3))) (let ((.def_17 (or A10 .def_16))) (let ((.def_18 (and .def_17 .def_15))) (let ((.def_19 (not .def_18))) (let ((.def_20 (not A11))) (let ((.def_21 (and A28 .def_20))) (let ((.def_22 (not A28))) (let ((.def_23 (not A26))) (let ((.def_24 (and .def_23 .def_22))) (let ((.def_25 (or .def_24 .def_21))) (let ((.def_26 (not .def_25))) (let ((.def_27 (or .def_26 .def_19))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and .def_28 .def_12))) (let ((.def_30 (not A29))) (let ((.def_31 (and .def_30 .def_16))) (let ((.def_32 (or A28 .def_22))) (let ((.def_33 (not .def_32))) (let ((.def_34 (or .def_33 .def_31))) (let ((.def_35 (not A9))) (let ((.def_36 (and A5 .def_35))) (let ((.def_37 (not A21))) (let ((.def_38 (or A12 .def_37))) (let ((.def_39 (not .def_38))) (let ((.def_40 (and .def_39 .def_36))) (let ((.def_41 (or .def_40 .def_34))) (let ((.def_42 (not .def_41))) (let ((.def_43 (or .def_35 A2))) (let ((.def_44 (not .def_43))) (let ((.def_45 (not A0))) (let ((.def_46 (or A21 .def_45))) (let ((.def_47 (not .def_46))) (let ((.def_48 (or .def_47 .def_44))) (let ((.def_49 (not .def_48))) (let ((.def_50 (= A5 A23))) (let ((.def_51 (and A9 A1))) (let ((.def_52 (not .def_51))) (let ((.def_53 (and .def_52 .def_50))) (let ((.def_54 (or .def_53 .def_49))) (let ((.def_55 (= .def_54 .def_42))) (let ((.def_56 (or .def_55 .def_29))) (let ((.def_57 (not A22))) (let ((.def_58 (= A12 .def_57))) (let ((.def_59 (and .def_35 .def_45))) (let ((.def_60 (and .def_59 .def_58))) (let ((.def_61 (not .def_60))) (let ((.def_62 (and .def_20 A2))) (let ((.def_63 (not A8))) (let ((.def_64 (and .def_63 A22))) (let ((.def_65 (or .def_64 .def_62))) (let ((.def_66 (or .def_65 .def_61))) (let ((.def_67 (not .def_66))) (let ((.def_68 (not A2))) (let ((.def_69 (or .def_30 .def_68))) (let ((.def_70 (and A27 .def_57))) (let ((.def_71 (not .def_70))) (let ((.def_72 (and .def_71 .def_69))) (let ((.def_73 (and A17 A24))) (let ((.def_74 (not .def_73))) (let ((.def_75 (and .def_35 A2))) (let ((.def_76 (not .def_75))) (let ((.def_77 (and .def_76 .def_74))) (let ((.def_78 (not .def_77))) (let ((.def_79 (= .def_78 .def_72))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or .def_80 .def_67))) (let ((.def_82 (not .def_81))) (let ((.def_83 (or A23 .def_68))) (let ((.def_84 (not .def_83))) (let ((.def_85 (and A1 A9))) (let ((.def_86 (not .def_85))) (let ((.def_87 (or .def_86 .def_84))) (let ((.def_88 (or A28 .def_1))) (let ((.def_89 (not .def_88))) (let ((.def_90 (not A17))) (let ((.def_91 (or .def_90 A12))) (let ((.def_92 (and .def_91 .def_89))) (let ((.def_93 (or .def_92 .def_87))) (let ((.def_94 (not .def_93))) (let ((.def_95 (not A18))) (let ((.def_96 (or .def_95 A15))) (let ((.def_97 (not A12))) (let ((.def_98 (and .def_30 .def_97))) (let ((.def_99 (not .def_98))) (let ((.def_100 (or .def_99 .def_96))) (let ((.def_101 (not .def_100))) (let ((.def_102 (and A1 .def_1))) (let ((.def_103 (and .def_90 A11))) (let ((.def_104 (and .def_103 .def_102))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_101))) (let ((.def_107 (not .def_106))) (let ((.def_108 (or .def_107 .def_94))) (let ((.def_109 (= .def_108 .def_82))) (let ((.def_110 (and .def_109 .def_56))) (let ((.def_111 (not .def_110))) .def_111)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)