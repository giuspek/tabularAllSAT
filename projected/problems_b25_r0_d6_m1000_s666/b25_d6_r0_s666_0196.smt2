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
(declare-fun A14 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A22))) (let ((.def_1 (not A5))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not A14))) (let ((.def_4 (and A19 .def_3))) (let ((.def_5 (not .def_4))) (let ((.def_6 (and .def_5 .def_2))) (let ((.def_7 (not .def_6))) (let ((.def_8 (not A17))) (let ((.def_9 (not A1))) (let ((.def_10 (and .def_9 .def_8))) (let ((.def_11 (and A18 A16))) (let ((.def_12 (= .def_11 .def_10))) (let ((.def_13 (not .def_12))) (let ((.def_14 (or .def_13 .def_7))) (let ((.def_15 (not .def_14))) (let ((.def_16 (not A21))) (let ((.def_17 (and .def_16 A5))) (let ((.def_18 (not .def_17))) (let ((.def_19 (or A5 .def_16))) (let ((.def_20 (not .def_19))) (let ((.def_21 (or .def_20 .def_18))) (let ((.def_22 (not A23))) (let ((.def_23 (or .def_22 A7))) (let ((.def_24 (not A20))) (let ((.def_25 (= .def_22 .def_24))) (let ((.def_26 (or .def_25 .def_23))) (let ((.def_27 (not .def_26))) (let ((.def_28 (or .def_27 .def_21))) (let ((.def_29 (and .def_28 .def_15))) (let ((.def_30 (not .def_29))) (let ((.def_31 (not A8))) (let ((.def_32 (or A22 .def_31))) (let ((.def_33 (not .def_32))) (let ((.def_34 (and .def_3 A1))) (let ((.def_35 (not .def_34))) (let ((.def_36 (and .def_35 .def_33))) (let ((.def_37 (or .def_24 A16))) (let ((.def_38 (and A5 A1))) (let ((.def_39 (not .def_38))) (let ((.def_40 (and .def_39 .def_37))) (let ((.def_41 (and .def_40 .def_36))) (let ((.def_42 (and A4 A24))) (let ((.def_43 (not .def_42))) (let ((.def_44 (not A4))) (let ((.def_45 (not A2))) (let ((.def_46 (or .def_45 .def_44))) (let ((.def_47 (or .def_46 .def_43))) (let ((.def_48 (not .def_47))) (let ((.def_49 (not A11))) (let ((.def_50 (or .def_49 A0))) (let ((.def_51 (not A19))) (let ((.def_52 (or .def_51 A22))) (let ((.def_53 (= .def_52 .def_50))) (let ((.def_54 (not .def_53))) (let ((.def_55 (or .def_54 .def_48))) (let ((.def_56 (or .def_55 .def_41))) (let ((.def_57 (and .def_56 .def_30))) (let ((.def_58 (not .def_57))) (let ((.def_59 (= .def_3 A21))) (let ((.def_60 (not .def_59))) (let ((.def_61 (and A16 A3))) (let ((.def_62 (not .def_61))) (let ((.def_63 (and .def_62 .def_60))) (let ((.def_64 (and A11 A11))) (let ((.def_65 (not .def_64))) (let ((.def_66 (not A6))) (let ((.def_67 (or .def_44 .def_66))) (let ((.def_68 (= .def_67 .def_65))) (let ((.def_69 (not .def_68))) (let ((.def_70 (and .def_69 .def_63))) (let ((.def_71 (and A3 A20))) (let ((.def_72 (not A24))) (let ((.def_73 (or .def_0 .def_72))) (let ((.def_74 (and .def_73 .def_71))) (let ((.def_75 (not .def_74))) (let ((.def_76 (and .def_31 A10))) (let ((.def_77 (or .def_31 A11))) (let ((.def_78 (or .def_77 .def_76))) (let ((.def_79 (not .def_78))) (let ((.def_80 (and .def_79 .def_75))) (let ((.def_81 (= .def_80 .def_70))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and .def_16 .def_72))) (let ((.def_84 (or .def_22 A9))) (let ((.def_85 (or .def_84 .def_83))) (let ((.def_86 (not .def_85))) (let ((.def_87 (and A0 A2))) (let ((.def_88 (or .def_3 .def_22))) (let ((.def_89 (and .def_88 .def_87))) (let ((.def_90 (or .def_89 .def_86))) (let ((.def_91 (or A11 .def_0))) (let ((.def_92 (and A16 .def_44))) (let ((.def_93 (not .def_92))) (let ((.def_94 (or .def_93 .def_91))) (let ((.def_95 (not A10))) (let ((.def_96 (or .def_95 .def_16))) (let ((.def_97 (not A9))) (let ((.def_98 (or .def_1 .def_97))) (let ((.def_99 (not .def_98))) (let ((.def_100 (or .def_99 .def_96))) (let ((.def_101 (not .def_100))) (let ((.def_102 (or .def_101 .def_94))) (let ((.def_103 (and .def_102 .def_90))) (let ((.def_104 (not .def_103))) (let ((.def_105 (or .def_104 .def_82))) (let ((.def_106 (or .def_105 .def_58))) (let ((.def_107 (not .def_106))) .def_107)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
