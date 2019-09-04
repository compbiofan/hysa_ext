use warnings;
use strict;

package graph;

sub rm_vertices{
    my ($g, $vertices) = @_;
    my $gg;
    foreach my $e (keys %$g){
        foreach my $e1 (keys %{$g->{$e}}){
            if(!defined $vertices->{$e} && !defined $vertices->{$e1}){
                $gg->{$e}->{$e1} = 1;
            }
        }
    }
    return $gg;
}

sub combine_graphs{
    my ($g1, $g2) = @_;
    my $g = $g1;
    foreach my $k1 (keys %$g2){
        foreach my $k2 (keys %{$g2->{$k1}}){
            $g->{$k1}->{$k2} = 1;
        }
    }
    return $g;
}
sub rm_edges{
    my ($g, $edges) = @_;
    my $gg;
    foreach my $e (keys %$g){
        foreach my $e1 (keys %{$g->{$e}}){
            if(!defined $edges->{$e}->{$e1} && !defined $edges->{$e1}->{$e}){
                $gg->{$e}->{$e1} = 1;
            }
        }
    }
    return $gg;
}
# output the number of nodes and edges
sub get_size{
    my ($g) = @_;
    my $h;
    my $edge_n = 0;
    foreach my $e (keys %$g){
        $h->{$e} = 1;
        foreach my $e1 (keys %{$g->{$e}}){
            $edge_n ++;
            $h->{$e1} = 1;
        }
    }
    my $node_n = scalar(keys %$h);
    return ($node_n, $edge_n);
}

sub print_m4{
    my ($g, $file) = @_;
    open fh_, ">$file" or die $!;
    foreach my $k (keys %$g){
        if($k !~ /\./){
            foreach my $kk (keys %{$g->{$k}}){
                print fh_ join("\t", $k, $kk) . "\n";
            }
        }
    }
    close fh_;
}

sub print_tree{
    my ($g) = @_;
    foreach my $k (keys %$g){
        foreach my $k1 (keys %{$g->{$k}}){
            print join("\t", $k, $k1) . "\n";
        }
    }
}

sub dfs_m{
    my ($g, $sources) = @_;
    # make the graph bidirectional
    foreach my $e1 (keys %$g){
        foreach my $e2 (keys %{$g->{$e1}}){
            $g->{$e2}->{$e1} = 1;
        }
    }
    my $debug = 0;
    my $id = 0;
    # ori -> new tree node coordination conversion
    my $index;
    # reverse
    my $index_rev;
    # the tree in the ordered new coordinate
    $index->{-1} = -1;
    $index_rev->{-1} = -1;
    my $gg;
    foreach my $source (keys %$sources){
        my $v;
        # the tree in the original coordinate
        my $ed;
        # for temporary record the edges in the order they are pushed in
        my @edd;
        my @s = ();
        push @s, $source;
        push @edd, "$source:-1";
        while(scalar(@s) != 0){
            my $c = pop @s;
            my $e = pop @edd;
            if(defined $v->{$c}){
                next;
            }
            # any node after this line should be recorded
            $v->{$c} = 1;
            $index->{$c}  = $id;
            $index_rev->{$id} = $c;
            my ($a, $b) = split(/:/, $e);
            $ed->{$a}->{$b} = 1;
            $id ++;
            foreach my $nb (keys %{$g->{$c}}){
                push @s, $nb;
                push @edd, "$nb:$c";
            }
        }
        if($debug == 1){
            print "Tree before conversion:\n";
            &print_tree($ed);
        }

        # return in the coordinate of dfs
        foreach my $n1 (keys %$ed){
            foreach my $n2 (keys %{$ed->{$n1}}){
                $gg->{$index->{$n1}}->{$index->{$n2}} = 1;
            }
        }
    }
    return ($gg, $index, $index_rev);
}
sub dfs{
    my ($g, $source) = @_;
    my $debug = 0;
    my @s = ();
    my $v;
    my $id = 0;
    # ori -> new tree node coordination conversion
    my $index;
    # reverse
    my $index_rev;
    # the tree in the original coordinate
    my $ed;
    # for temporary record the edges in the order they are pushed in
    my @edd;
    # the tree in the ordered new coordinate
    my $gg;
    push @s, $source;
    push @edd, "$source:-1";
    $index->{-1} = -1;
    $index_rev->{-1} = -1;
    while(scalar(@s) != 0){
        my $c = pop @s;
        my $e = pop @edd;
        if(defined $v->{$c}){
            next;
        }
        # any node after this line should be recorded
        $v->{$c} = 1;
        $index->{$c}  = $id;
        $index_rev->{$id} = $c;
        my ($a, $b) = split(/:/, $e);
        $ed->{$a}->{$b} = 1;
        $id ++;
        foreach my $nb (keys %{$g->{$c}}){
            push @s, $nb;
            push @edd, "$nb:$c";
        }
    }
    if($debug == 1){
        print "Tree before conversion:\n";
        &print_tree($ed);
    }

    # return in the coordinate of dfs
    foreach my $n1 (keys %$ed){
        foreach my $n2 (keys %{$ed->{$n1}}){
            $gg->{$index->{$n1}}->{$index->{$n2}} = 1;
        }
    }
    return ($gg, $index, $index_rev);
}

sub find_missing_edge{
    my ($g, $gg, $index) =  @_;
    my $e;

    foreach my $c (keys %$g){
        # get gg's coordinate of c in g
        if(!defined $index->{$c}){
            my $a = 1;
        }
        my $c_ = $index->{$c};
        foreach my $n (keys %{$g->{$c}}){
            my $n_ = $index->{$n};
            if(defined $gg->{$c_}->{$n_} || defined $gg->{$n_}->{$c_}){
                next;
            }
            if($c_ > $n_){
                $e->{$n_}->{$c_} = 1;
            }
            else{
                $e->{$c_}->{$n_} = 1;
            }
        }
    }
    return $e;
}

1;


