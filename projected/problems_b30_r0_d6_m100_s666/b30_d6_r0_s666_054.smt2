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
(declare-fun A25 () Bool)
(declare-fun A27 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A6))) (let ((.def_1 (not A4))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not .def_2))) (let ((.def_4 (not A7))) (let ((.def_5 (or A19 .def_4))) (let ((.def_6 (not .def_5))) (let ((.def_7 (or .def_6 .def_3))) (let ((.def_8 (not A5))) (let ((.def_9 (or A8 .def_8))) (let ((.def_10 (not .def_9))) (let ((.def_11 (or A2 A23))) (let ((.def_12 (or .def_11 .def_10))) (let ((.def_13 (or .def_12 .def_7))) (let ((.def_14 (not A20))) (let ((.def_15 (and .def_14 .def_0))) (let ((.def_16 (and A0 A21))) (let ((.def_17 (and .def_16 .def_15))) (let ((.def_18 (not .def_17))) (let ((.def_19 (not A12))) (let ((.def_20 (or A24 .def_19))) (let ((.def_21 (not .def_20))) (let ((.def_22 (or A11 A15))) (let ((.def_23 (or .def_22 .def_21))) (let ((.def_24 (not .def_23))) (let ((.def_25 (or .def_24 .def_18))) (let ((.def_26 (and .def_25 .def_13))) (let ((.def_27 (and A3 A7))) (let ((.def_28 (not .def_27))) (let ((.def_29 (not A2))) (let ((.def_30 (or A20 .def_29))) (let ((.def_31 (and .def_30 .def_28))) (let ((.def_32 (not .def_31))) (let ((.def_33 (and A19 A18))) (let ((.def_34 (or A24 .def_4))) (let ((.def_35 (or .def_34 .def_33))) (let ((.def_36 (or .def_35 .def_32))) (let ((.def_37 (and A4 A21))) (let ((.def_38 (not .def_37))) (let ((.def_39 (not A15))) (let ((.def_40 (not A19))) (let ((.def_41 (and .def_40 .def_39))) (let ((.def_42 (not .def_41))) (let ((.def_43 (= .def_42 .def_38))) (let ((.def_44 (and .def_14 .def_40))) (let ((.def_45 (not .def_44))) (let ((.def_46 (= .def_14 A12))) (let ((.def_47 (not .def_46))) (let ((.def_48 (or .def_47 .def_45))) (let ((.def_49 (not .def_48))) (let ((.def_50 (or .def_49 .def_43))) (let ((.def_51 (not .def_50))) (let ((.def_52 (= .def_51 .def_36))) (let ((.def_53 (or .def_52 .def_26))) (let ((.def_54 (not .def_53))) (let ((.def_55 (not A25))) (let ((.def_56 (and .def_55 A11))) (let ((.def_57 (not .def_56))) (let ((.def_58 (or .def_0 A29))) (let ((.def_59 (not .def_58))) (let ((.def_60 (or .def_59 .def_57))) (let ((.def_61 (not .def_60))) (let ((.def_62 (and .def_29 A14))) (let ((.def_63 (and A17 A6))) (let ((.def_64 (not .def_63))) (let ((.def_65 (and .def_64 .def_62))) (let ((.def_66 (or .def_65 .def_61))) (let ((.def_67 (not .def_66))) (let ((.def_68 (or A4 A27))) (let ((.def_69 (or A0 A23))) (let ((.def_70 (not .def_69))) (let ((.def_71 (or .def_70 .def_68))) (let ((.def_72 (not .def_71))) (let ((.def_73 (and .def_8 A4))) (let ((.def_74 (not .def_73))) (let ((.def_75 (not A10))) (let ((.def_76 (and .def_1 .def_75))) (let ((.def_77 (= .def_76 .def_74))) (let ((.def_78 (not .def_77))) (let ((.def_79 (= .def_78 .def_72))) (let ((.def_80 (not .def_79))) (let ((.def_81 (or .def_80 .def_67))) (let ((.def_82 (not A22))) (let ((.def_83 (or A25 .def_82))) (let ((.def_84 (not .def_83))) (let ((.def_85 (not A16))) (let ((.def_86 (or A8 .def_85))) (let ((.def_87 (or .def_86 .def_84))) (let ((.def_88 (and .def_39 .def_75))) (let ((.def_89 (not .def_88))) (let ((.def_90 (or A3 A9))) (let ((.def_91 (not .def_90))) (let ((.def_92 (and .def_91 .def_89))) (let ((.def_93 (or .def_92 .def_87))) (let ((.def_94 (and .def_82 A2))) (let ((.def_95 (or .def_55 .def_8))) (let ((.def_96 (and .def_95 .def_94))) (let ((.def_97 (not .def_96))) (let ((.def_98 (not A24))) (let ((.def_99 (not A27))) (let ((.def_100 (and .def_99 .def_98))) (let ((.def_101 (not .def_100))) (let ((.def_102 (or A0 A8))) (let ((.def_103 (not .def_102))) (let ((.def_104 (and .def_103 .def_101))) (let ((.def_105 (and .def_104 .def_97))) (let ((.def_106 (not .def_105))) (let ((.def_107 (and .def_106 .def_93))) (let ((.def_108 (not .def_107))) (let ((.def_109 (and .def_108 .def_81))) (let ((.def_110 (not .def_109))) (let ((.def_111 (and .def_110 .def_54))) .def_111)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
