(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A24 () Bool)
(declare-fun A25 () Bool)
(declare-fun A26 () Bool)
(declare-fun A27 () Bool)
(declare-fun A28 () Bool)
(declare-fun A29 () Bool)
(assert (let ((.def_0 (not A24))) (let ((.def_1 (or A4 .def_0))) (let ((.def_2 (not .def_1))) (let ((.def_3 (or A28 A17))) (let ((.def_4 (not .def_3))) (let ((.def_5 (or .def_4 .def_2))) (let ((.def_6 (not A29))) (let ((.def_7 (not A13))) (let ((.def_8 (or .def_7 .def_6))) (let ((.def_9 (not A26))) (let ((.def_10 (or .def_9 A27))) (let ((.def_11 (or .def_10 .def_8))) (let ((.def_12 (not .def_11))) (let ((.def_13 (or .def_12 .def_5))) (let ((.def_14 (not A16))) (let ((.def_15 (or A20 .def_14))) (let ((.def_16 (not A11))) (let ((.def_17 (and .def_16 A6))) (let ((.def_18 (not .def_17))) (let ((.def_19 (= .def_18 .def_15))) (let ((.def_20 (not A14))) (let ((.def_21 (and .def_0 .def_20))) (let ((.def_22 (not .def_21))) (let ((.def_23 (not A4))) (let ((.def_24 (or A3 .def_23))) (let ((.def_25 (or .def_24 .def_22))) (let ((.def_26 (or .def_25 .def_19))) (let ((.def_27 (or .def_26 .def_13))) (let ((.def_28 (not .def_27))) (let ((.def_29 (not A9))) (let ((.def_30 (not A6))) (let ((.def_31 (and .def_30 .def_29))) (let ((.def_32 (not .def_31))) (let ((.def_33 (and A0 A16))) (let ((.def_34 (not .def_33))) (let ((.def_35 (and .def_34 .def_32))) (let ((.def_36 (not A25))) (let ((.def_37 (or .def_36 .def_30))) (let ((.def_38 (not .def_37))) (let ((.def_39 (not A1))) (let ((.def_40 (and .def_39 A28))) (let ((.def_41 (and .def_40 .def_38))) (let ((.def_42 (or .def_41 .def_35))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A10))) (let ((.def_45 (and A22 .def_44))) (let ((.def_46 (not .def_45))) (let ((.def_47 (and A17 A28))) (let ((.def_48 (or .def_47 .def_46))) (let ((.def_49 (and A5 A9))) (let ((.def_50 (not .def_49))) (let ((.def_51 (not A7))) (let ((.def_52 (and .def_51 .def_0))) (let ((.def_53 (not .def_52))) (let ((.def_54 (= .def_53 .def_50))) (let ((.def_55 (not .def_54))) (let ((.def_56 (and .def_55 .def_48))) (let ((.def_57 (not .def_56))) (let ((.def_58 (and .def_57 .def_43))) (let ((.def_59 (or .def_58 .def_28))) (let ((.def_60 (and A22 A27))) (let ((.def_61 (not .def_60))) (let ((.def_62 (not A28))) (let ((.def_63 (and .def_14 .def_62))) (let ((.def_64 (and .def_63 .def_61))) (let ((.def_65 (not A18))) (let ((.def_66 (= A10 .def_65))) (let ((.def_67 (not A0))) (let ((.def_68 (or A21 .def_67))) (let ((.def_69 (or .def_68 .def_66))) (let ((.def_70 (not .def_69))) (let ((.def_71 (and .def_70 .def_64))) (let ((.def_72 (not .def_71))) (let ((.def_73 (= A17 .def_51))) (let ((.def_74 (and .def_23 .def_9))) (let ((.def_75 (and .def_74 .def_73))) (let ((.def_76 (or .def_65 A11))) (let ((.def_77 (not .def_76))) (let ((.def_78 (or A3 .def_51))) (let ((.def_79 (or .def_78 .def_77))) (let ((.def_80 (not .def_79))) (let ((.def_81 (and .def_80 .def_75))) (let ((.def_82 (or .def_81 .def_72))) (let ((.def_83 (not .def_82))) (let ((.def_84 (or A25 .def_39))) (let ((.def_85 (not A5))) (let ((.def_86 (or .def_85 A9))) (let ((.def_87 (and .def_86 .def_84))) (let ((.def_88 (or .def_44 A19))) (let ((.def_89 (not A22))) (let ((.def_90 (or .def_89 .def_67))) (let ((.def_91 (and .def_90 .def_88))) (let ((.def_92 (and .def_91 .def_87))) (let ((.def_93 (not .def_92))) (let ((.def_94 (or A29 A6))) (let ((.def_95 (and .def_51 A6))) (let ((.def_96 (and .def_95 .def_94))) (let ((.def_97 (not .def_96))) (let ((.def_98 (not A20))) (let ((.def_99 (= .def_98 A14))) (let ((.def_100 (or A18 .def_6))) (let ((.def_101 (and .def_100 .def_99))) (let ((.def_102 (and .def_101 .def_97))) (let ((.def_103 (not .def_102))) (let ((.def_104 (or .def_103 .def_93))) (let ((.def_105 (not .def_104))) (let ((.def_106 (and .def_105 .def_83))) (let ((.def_107 (or .def_106 .def_59))) .def_107)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
