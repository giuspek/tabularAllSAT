(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
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
(assert (let ((.def_0 (not A5))) (let ((.def_1 (or .def_0 A20))) (let ((.def_2 (not .def_1))) (let ((.def_3 (or A14 A10))) (let ((.def_4 (not .def_3))) (let ((.def_5 (or .def_4 .def_2))) (let ((.def_6 (not A2))) (let ((.def_7 (= A26 .def_6))) (let ((.def_8 (not .def_7))) (let ((.def_9 (not A23))) (let ((.def_10 (not A6))) (let ((.def_11 (= .def_10 .def_9))) (let ((.def_12 (not .def_11))) (let ((.def_13 (and .def_12 .def_8))) (let ((.def_14 (and .def_13 .def_5))) (let ((.def_15 (not .def_14))) (let ((.def_16 (or A2 A12))) (let ((.def_17 (not .def_16))) (let ((.def_18 (not A15))) (let ((.def_19 (= .def_18 A28))) (let ((.def_20 (not .def_19))) (let ((.def_21 (= .def_20 .def_17))) (let ((.def_22 (not .def_21))) (let ((.def_23 (not A13))) (let ((.def_24 (or A5 .def_23))) (let ((.def_25 (and .def_10 A21))) (let ((.def_26 (not .def_25))) (let ((.def_27 (and .def_26 .def_24))) (let ((.def_28 (not .def_27))) (let ((.def_29 (and .def_28 .def_22))) (let ((.def_30 (or .def_29 .def_15))) (let ((.def_31 (not .def_30))) (let ((.def_32 (or A22 A18))) (let ((.def_33 (or A11 A29))) (let ((.def_34 (not .def_33))) (let ((.def_35 (and .def_34 .def_32))) (let ((.def_36 (or A10 A16))) (let ((.def_37 (and A17 A27))) (let ((.def_38 (and .def_37 .def_36))) (let ((.def_39 (not .def_38))) (let ((.def_40 (= .def_39 .def_35))) (let ((.def_41 (not .def_40))) (let ((.def_42 (= A6 A1))) (let ((.def_43 (not .def_42))) (let ((.def_44 (= .def_18 A17))) (let ((.def_45 (not .def_44))) (let ((.def_46 (or .def_45 .def_43))) (let ((.def_47 (not .def_46))) (let ((.def_48 (not A16))) (let ((.def_49 (and A24 .def_48))) (let ((.def_50 (not A28))) (let ((.def_51 (and .def_6 .def_50))) (let ((.def_52 (= .def_51 .def_49))) (let ((.def_53 (not .def_52))) (let ((.def_54 (and .def_53 .def_47))) (let ((.def_55 (not .def_54))) (let ((.def_56 (= .def_55 .def_41))) (let ((.def_57 (and .def_56 .def_31))) (let ((.def_58 (not A26))) (let ((.def_59 (not A19))) (let ((.def_60 (or .def_59 .def_58))) (let ((.def_61 (or .def_36 .def_60))) (let ((.def_62 (not .def_61))) (let ((.def_63 (not A14))) (let ((.def_64 (not A0))) (let ((.def_65 (or .def_64 .def_63))) (let ((.def_66 (not .def_65))) (let ((.def_67 (and A3 .def_6))) (let ((.def_68 (and .def_67 .def_66))) (let ((.def_69 (not .def_68))) (let ((.def_70 (and .def_69 .def_62))) (let ((.def_71 (not .def_70))) (let ((.def_72 (not A9))) (let ((.def_73 (and .def_72 A6))) (let ((.def_74 (not .def_73))) (let ((.def_75 (and A25 A29))) (let ((.def_76 (and .def_75 .def_74))) (let ((.def_77 (or A11 .def_18))) (let ((.def_78 (not .def_77))) (let ((.def_79 (and A26 A20))) (let ((.def_80 (or .def_79 .def_78))) (let ((.def_81 (and .def_80 .def_76))) (let ((.def_82 (or .def_81 .def_71))) (let ((.def_83 (not .def_82))) (let ((.def_84 (and A18 .def_18))) (let ((.def_85 (or A5 A23))) (let ((.def_86 (not .def_85))) (let ((.def_87 (= .def_86 .def_84))) (let ((.def_88 (not .def_87))) (let ((.def_89 (not A24))) (let ((.def_90 (and .def_89 .def_48))) (let ((.def_91 (and A28 .def_9))) (let ((.def_92 (not .def_91))) (let ((.def_93 (or .def_92 .def_90))) (let ((.def_94 (and .def_93 .def_88))) (let ((.def_95 (not .def_94))) (let ((.def_96 (and A22 .def_48))) (let ((.def_97 (not .def_96))) (let ((.def_98 (or A12 A28))) (let ((.def_99 (or .def_98 .def_97))) (let ((.def_100 (not .def_99))) (let ((.def_101 (not A3))) (let ((.def_102 (or .def_64 .def_101))) (let ((.def_103 (not A8))) (let ((.def_104 (or .def_103 .def_101))) (let ((.def_105 (not .def_104))) (let ((.def_106 (= .def_105 .def_102))) (let ((.def_107 (= .def_106 .def_100))) (let ((.def_108 (and .def_107 .def_95))) (let ((.def_109 (not .def_108))) (let ((.def_110 (and .def_109 .def_83))) (let ((.def_111 (not .def_110))) (let ((.def_112 (or .def_111 .def_57))) (let ((.def_113 (not .def_112))) .def_113)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)